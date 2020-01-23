#!/usr/bin/env python3
import itertools
import logging
from functools import reduce
logging.basicConfig(level=logging.DEBUG)
import os
import re
import shutil
from io import StringIO
from argparse import ArgumentParser
from copy import deepcopy
from collections import namedtuple, deque, defaultdict
from subprocess import check_output
from tempfile import NamedTemporaryFile
from textwrap import dedent
from urllib.parse import urlparse

import newick
import pytest
from attr import attrs, attrib
from sonLib.bioio import reverseComplement, fastaRead
from toil.common import Toil
from toil.job import Job

logger = logging.getLogger(__name__)

def get_ascii_tmp_file(delete=True):
    "create a temporary file with ASCII encode"
    return NamedTemporaryFile(mode="w", encoding='ascii', delete=delete)

class LiftedContextManager(object):
    """Transforms an iterable of context managers into a list of opened managers."""
    def __init__(self, iterable):
        self.iterable = iterable

    def __enter__(self):
        return [i.__enter__() for i in self.iterable]

    def __exit__(self, exc_type, exc_value, traceback):
        retval = False
        for contextmanager in self.iterable:
            retval = contextmanager.__exit__(exc_type, exc_value, traceback) or retval
        return retval

def call(cmd):
    """
    Call the given command and return the output.

    Why not just use check_output directly? This is in case we need to
    use Docker/Singularity at some point.
    """
    return check_output(cmd)

def get_hal_tree(hal_path):
    """Get the species tree of the given hal file."""
    return call(["halStats", "--tree", hal_path])

def find_node_by_name(tree, name):
    """Get a single node in the tree with the given name, or None if no such node exists."""
    for node in tree.walk():
        if node.name == name:
            return node
    return None

def clone_node(node):
    """Get a new Node with the same properties, but no ancestor or children."""
    node = newick.Node(node.name, str(node.length))
    return node

def reroot_tree(tree, new_root_name):
    """
    Reroot the tree at the node with name new_root_name.

    The root node is kept even if it becomes superfluous.
    """
    node = find_node_by_name(tree, new_root_name)
    if node is None:
        raise ValueError("Node {} not found in tree".format(new_root_name))
    if node.ancestor is None:
        # Already the root
        return deepcopy(node)

    def clone_flipped_tree(cur_node, old_node, node_to_add_to, branch_length):
        cloned_node = clone_node(cur_node)
        node_to_add_to.add_descendant(cloned_node)
        cloned_node.length = branch_length
        for child in cur_node.descendants:
            if child != old_node:
                cloned_child = deepcopy(child)
                cloned_node.add_descendant(cloned_child)
        if cur_node.ancestor is not None:
            clone_flipped_tree(cur_node.ancestor, cur_node, cloned_node, cur_node.length)

    new_root = deepcopy(node)
    new_root.length = None
    new_root.ancestor = None
    clone_flipped_tree(node.ancestor, node, new_root, node.length)
    return new_root

def prune_tree(tree, leaves_to_keep):
    """
    Remove any nodes except the given leaves and their
    ancestors. Ancestral nodes are kept even if they only have a single
    child.
    """
    nodes = tree.walk(mode='postorder')
    for node in nodes:
        if node.is_leaf and node.name not in leaves_to_keep:
            node.ancestor.descendants.remove(node)

def start_job(job, hal_id, refGenome, opts):
    """Set up the structure of the pipeline."""
    hal = hal_id
    # Newick representation of the HAL species tree.
    newick_string = get_hal_tree(hal)
    job.fileStore.logToMaster("Newick string: %s" % (newick_string))
    tree = newick.loads(newick_string)[0]
    rerooted = reroot_tree(tree, refGenome)
    job.fileStore.logToMaster("Rerooted newick string: %s" % (newick.dumps([rerooted])))
    if opts.targetGenomes is not None:
        # We don't need the alignment to all genomes, just a subset.
        prune_tree(rerooted, opts.targetGenomes)
        job.fileStore.logToMaster("Pruned newick string: %s" % newick.dumps(rerooted))

    def setup_jobs(node):
        """Recursively set up jobs for this node and its children."""
        prev_data = [setup_jobs(child) for child in node.descendants]
        # At this point all of the jobs for the lower parts of the tree have been set up.
        lifted_data = [prev_lifted for _, prev_lifted in prev_data]
        merge_job = job.wrapJobFn(merge_blocks_job, node.name, [n.name for n in node.descendants], lifted_data, hal_id, opts)
        for prev_job, _ in prev_data:
            prev_job.addFollowOn(merge_job)
        if node.is_leaf:
            job.addChild(merge_job)
        if node.ancestor is None:
            return merge_job.rv()
        else:
            # Find out whether we have to lift up or down
            original_node = find_node_by_name(tree, node.name)
            if original_node.ancestor is None or node.ancestor.name != original_node.ancestor.name:
                lift_down_job = merge_job.addFollowOnJobFn(lift_job, 'down', node.name, node.ancestor.name, merge_job.rv(), hal_id, opts)
                return lift_down_job, lift_down_job.rv()
            else:
                lift_up_job = merge_job.addFollowOnJobFn(lift_job, 'up', node.name, node.ancestor.name, merge_job.rv(), hal_id, opts)
                return lift_up_job, lift_up_job.rv()

    blocks_on_ref = setup_jobs(rerooted)

    all_genomes = [node.name for node in tree.walk()]

    return job.addFollowOnJobFn(maf_export_job, hal, all_genomes, blocks_on_ref, opts).rv()

def split_blocks_by_chrom(blocks_file, get_temp_path):
    """
    Split a blocks file into several files, one per chromosome.

    get_temp_path is a function that returns a new temporary file path.

    Returns a dictionary, {chrom_name: blockfile_path}.
    """
    paths = {}
    cur_chrom = None
    output_path = get_temp_path()
    output = open(output_path, 'w')
    while True:
        block = Block.read_next_from_file(blocks_file)
        if block is None:
            break
        if cur_chrom is None:
            cur_chrom = block.chrom
            paths[cur_chrom] = output_path
        elif block.chrom != cur_chrom:
            # Done with this chromosome
            output.close()
            output_path = get_temp_path()
            cur_chrom = block.chrom
            assert cur_chrom not in paths
            paths[cur_chrom] = output_path
            output = open(output_path, 'w')
        output.write(str(block))
    return paths

def split_alignments_by_chrom(alignments_file, get_temp_path):
    """Split an alignments file into several files, one per chromosome.

    get_temp_path is a function that returns a new temporary file path.

    Returns a dictionary, {chrom_name: alignmentfile_path}."""
    paths = {}
    cur_chrom = None
    output_path = get_temp_path()
    output = open(output_path, 'w')
    for line in alignments_file:
        fields = line.split()
        line_chrom = fields[0]
        if line_chrom != cur_chrom:
            output.close()
            output_path = get_temp_path()
            cur_chrom = line_chrom
            paths[cur_chrom] = output_path
            output = open(output_path, 'w')
        output.write(line)
    return paths

def group_by_size(files_by_chroms, chrom_sizes, get_temp_path, target_size=50000000):
    """
    Group a series of files (each indexed by chromosome) into larger series of files
    with total sequence approximately target_size.

    Returns two lists like [chrom_sizes_1, chrom_sizes_2], [[alignments_file_1, blocks_file_1], ...].
    """
    def merge(files_by_chroms, chroms_to_merge):
        # Ensure the output will be generated in sorted order
        chroms_to_merge.sort()

        paths = [get_temp_path() for _ in range(len(files_by_chroms))]
        with LiftedContextManager([open(p, 'w') for p in paths]) as outputs:
            for chrom in chroms_to_merge:
                for files_by_chrom, output in zip(files_by_chroms, outputs):
                    if chrom not in files_by_chrom:
                        # No file corresponding to this chrom
                        continue
                    with open(files_by_chrom[chrom]) as input:
                        shutil.copyfileobj(input, output)
        return paths

    grouped_pathss = []
    chrom_sizes_subsets = []
    chromsAndSizes = sorted(list(chrom_sizes.items()), key=lambda cs: cs[1])
    running_size = 0
    chroms_to_merge = []
    for chrom, size in chromsAndSizes:
        if running_size > 0 and running_size + size > target_size:
            grouped_paths = merge(files_by_chroms, chroms_to_merge)
            grouped_pathss.append(grouped_paths)
            chrom_sizes_subsets.append({c: chrom_sizes[c] for c in chroms_to_merge})
            chroms_to_merge = []
            running_size = 0
        chroms_to_merge.append(chrom)
        running_size += size
    grouped_paths = merge(files_by_chroms, chroms_to_merge)
    chrom_sizes_subsets.append({c: chrom_sizes[c] for c in chroms_to_merge})
    grouped_pathss.append(grouped_paths)

    assert len(chrom_sizes_subsets) == len(grouped_pathss)
    return chrom_sizes_subsets, grouped_pathss

def lift_job(job, up_or_down, genome, other_genome, blocks_id, hal_id, opts):
    """Lift this genome's blocks file to another genome."""
    job.fileStore.logToMaster("Lifting from {} to {}".format(genome, other_genome))

    # Get the alignments / blocks, separated into different files by chromosome.
    hal = hal_id
    alignment_path = job.fileStore.getLocalTempFile()
    if up_or_down == 'up':
        get_alignment_to_parent(hal, genome, alignment_path)
    else:
        assert up_or_down == 'down'
        get_alignment_to_child(hal, other_genome, alignment_path)
    with open(alignment_path) as alignments:
        alignments_by_chrom = split_alignments_by_chrom(alignments, job.fileStore.getLocalTempFile)
    blocks_path = job.fileStore.readGlobalFile(blocks_id)
    with open(blocks_path) as blocks:
        blocks_by_chrom = split_blocks_by_chrom(blocks, job.fileStore.getLocalTempFile)

    # Group up the files into larger chunks, each of which will be handled by a separate job.
    chrom_sizes = get_chrom_sizes(hal, genome)
    chrom_sizes_subsets, grouped_pathss = group_by_size([alignments_by_chrom, blocks_by_chrom], chrom_sizes, job.fileStore.getLocalTempFile)

    # Set up the jobs to handle each chunk.
    subset_jobs = []
    for chrom_sizes_subset, grouped_paths in zip(chrom_sizes_subsets, grouped_pathss):
        alignments_path = grouped_paths[0]
        blocks_path = grouped_paths[1]
        alignment_id = job.fileStore.writeGlobalFile(alignments_path)
        blocks_id = job.fileStore.writeGlobalFile(blocks_path)
        subset_job = job.addChildJobFn(lift_subset_job, alignment_id, blocks_id, other_genome)
        subset_jobs.append(subset_job)

    # Set up the jobs to handle the output of the chunked jobs (concatenating and sorting).
    cat_job = job.addChildJobFn(concatenate_job, [j.rv() for j in subset_jobs])
    for subset_job in subset_jobs:
        subset_job.addChild(cat_job)
    output_id = cat_job.addChildJobFn(sort_blocks_job, cat_job.rv()).rv()
    return output_id

def lift_subset_job(job, alignment_id, blocks_id, parent_genome):
    alignment_path = alignment_id
    blocks_path = job.fileStore.readGlobalFile(blocks_id)
    lifted_path = job.fileStore.getLocalTempFile()
    with open(alignment_path) as alignment, open(blocks_path) as blocks, open(lifted_path, 'w') as lifted:
        lift_blocks(alignment, blocks, parent_genome, lifted)
    return job.fileStore.writeGlobalFile(lifted_path)

def sort_blocks_job(job, blocks_id):
    blocks_path = job.fileStore.readGlobalFile(blocks_id)
    sorted_path = job.fileStore.getLocalTempFile()
    with open(blocks_path) as blocks, open(sorted_path, 'w') as sorted:
        sort_blocks(blocks, sorted)
    return job.fileStore.writeGlobalFile(sorted_path)

def merge_blocks_job(job, genome, child_names, block_ids, hal_id, opts):
    """Combine the lifted blocks and add in the current genome's sequence."""
    job.fileStore.logToMaster("Merging blocks on {} (from child genomes {})".format(genome, child_names))
    hal = hal_id
    if len(block_ids) == 0:
        # Leaf genome.
        assert len(child_names) == len(block_ids) == 0
        output_path = job.fileStore.getLocalTempFile()
        get_leaf_blocks(hal, genome, output_path)
        output_id = job.fileStore.writeGlobalFile(output_path)
    else:
        # Have several block_ids to merge.
        assert len(child_names) == len(block_ids)
        assert len(child_names) > 0

        chrom_sizes = get_chrom_sizes(hal, genome)

        block_paths = list(map(job.fileStore.readGlobalFile, block_ids))
        split_blocks = [split_blocks_by_chrom(open(block_path), job.fileStore.getLocalTempFile) for block_path in block_paths]
        chrom_sizes_subsets, grouped_pathss = group_by_size(split_blocks, chrom_sizes, job.fileStore.getLocalTempFile)
        subset_jobs = []
        for chrom_sizes_subset, grouped_paths in zip(chrom_sizes_subsets, grouped_pathss):
            grouped_ids = [job.fileStore.writeGlobalFile(grouped_path) for grouped_path in grouped_paths]
            subset_job = job.addChildJobFn(merge_blocks_merge_subset_job, genome, child_names, grouped_ids, chrom_sizes_subset, opts)
            subset_jobs.append(subset_job)

        cat_job = job.addChildJobFn(concatenate_job, [j.rv() for j in subset_jobs])
        for subset_job in subset_jobs:
            subset_job.addChild(cat_job)
        output_id = cat_job.addChildJobFn(merge_blocks_maximize_block_length_job, genome, cat_job.rv(), hal_id, opts).rv()
    return output_id

def merge_blocks_merge_subset_job(job, genome, child_names, block_ids, chrom_sizes, opts):
    """Merges a small subset of the lifted alignments from children/parent."""
    block_paths = [job.fileStore.readGlobalFile(block_id) for block_id in block_ids]
    merged_path = job.fileStore.getLocalTempFile()
    with LiftedContextManager(list(map(open, block_paths))) as blocks, open(merged_path, 'w') as output:
        merge_child_blocks(genome, chrom_sizes, child_names, blocks, output)
    output_id = job.fileStore.writeGlobalFile(merged_path)
    if opts.intermediateResultsUrl is not None:
        job.fileStore.exportFile(output_id, opts.intermediateResultsUrl + "{}-merged-subset-{}.blocks".format(genome, list(chrom_sizes.keys())[0]))
    return output_id

def merge_blocks_maximize_block_length_job(job, genome, merged_blocks_id, hal_id, opts):
    """Maximizes the block length of a set of merged blocks."""
    hal = hal_id
    ref_sequence = get_sequence(hal, genome)
    merged_path = job.fileStore.readGlobalFile(merged_blocks_id)
    # Sort the blocks (they are probably in a somewhat random order).
    sorted_path = job.fileStore.getLocalTempFile()
    with open(merged_path) as merged_blocks, open(sorted_path, 'w') as sorted_blocks:
        sort_blocks(merged_blocks, sorted_blocks)

    output_path = job.fileStore.getLocalTempFile()
    with open(sorted_path) as sorted_blocks, open(output_path, 'w') as output:
        maximize_block_length(sorted_blocks, ref_sequence, output)
    output_id = job.fileStore.writeGlobalFile(output_path)
    if opts.intermediateResultsUrl is not None:
        job.fileStore.exportFile(output_id, opts.intermediateResultsUrl + "{}-merged.blocks".format(genome))
    return output_id

def concatenate_job(job, input_ids):
    """This job just concatenates the input files together and returns them."""
    out_path = job.fileStore.getLocalTempFile()
    with open(out_path, 'w') as output:
        for input_id in input_ids:
            input_path = job.fileStore.readGlobalFile(input_id)
            with open(input_path) as input:
                shutil.copyfileobj(input, output)
    return job.fileStore.writeGlobalFile(out_path)

def find_split_point(blocks):
    """
    Find the position at which blocks stop sharing overlap.
    """
    return min([(b.chrom, b.end) for b in blocks if b is not None])[1]

def merged_dup_stream(file, read_next):
    """
    Stream blocks/alignments from a file, merging overlapping dups.

    Given a sorted file containing block/alignment objects (retrieved
    via the read_next function), this will eagerly split at any points
    of overlap and merge objects that belong to the same interval.
    """
    # NB: we refer to "blocks" here, but since AlignmentGroup responds to the
    # same calls as Block, we can use this function on alignments as well.
    queued_blocks = deque()
    while True:
        # Get our next block to consider.
        if len(queued_blocks) > 0:
            block = queued_blocks.popleft()
        else:
            block = read_next(file)
        if block is None:
            # End of file.
            break
        # Check for overlapping blocks in our queue.
        dup_blocks = []
        while len(queued_blocks) > 0:
            if block.overlaps(queued_blocks[0]) and all(b.overlaps(queued_blocks[0]) for b in dup_blocks):
                dup_blocks.append(queued_blocks.popleft())
            else:
                break
        # Check for overlapping blocks in our file.
        while True:
            next_block = read_next(file)
            if next_block is None:
                break
            if block.overlaps(next_block) and all(b.overlaps(next_block) for b in dup_blocks):
                dup_blocks.append(next_block)
            else:
                queued_blocks.append(next_block)
                break

        # Go through the overlapping blocks (if any) and perform the actual
        # merging / splitting.
        if len(dup_blocks) > 0:
            if not all(b.start == block.start for b in dup_blocks):
                # Not every block starts at the same point, so we look for where
                # the overlap *starts* and split there (so we can output the
                # part of the block that doesn't share overlap).
                split_point = min([b.start for b in dup_blocks if b.start != block.start])
            else:
                # In this case, every block starts at the same point, so we need
                # to look for where the overlap *stops*.
                split_point = min([b.end for b in [block] + dup_blocks])
            leftovers = []
            merged_block = None
            for dup_block in [block] + dup_blocks:
                if dup_block.start < split_point < dup_block.end:
                    # We have a block that overlaps the splitting point. Split
                    # it and save the rest for a later iteration.
                    dup_block, leftover = dup_block.split(split_point - dup_block.start)
                    leftovers.append(leftover)
                if merged_block is None:
                    # We share overlap with the emitted block by definition
                    # (this is always some part of the initial block).
                    merged_block = dup_block
                elif dup_block.start < split_point:
                    # We already have a block to emit, so we have to add in the
                    # information from this block.
                    merged_block = merged_block.merge(dup_block)
                else:
                    # No overlap with the emitted block. Save this for later.
                    leftovers.append(dup_block)
            block = merged_block
            # Save the blocks that we didn't use for later, maintaining sorted order.
            queued_blocks.extendleft(reversed(leftovers))
        logger.debug('yielding block %s from dup stream', block)
        yield block

def merge_child_blocks(genome, chrom_sizes, child_names, block_files, output):
    """Merge several block files on the same reference into one file.

    All input block files must be sorted. This merges the alignment information
    from the blocks "vertically" so that there is exactly 1 block per reference
    position, adding in a single-degree block if no input blocks covered a
    reference position.

    The reference sequence is left as 'X', to be filled in at a later time.

    """
    def get_smallest_block(cur_blocks):
        return reduce(lambda a, i: a if a.first < i.first else i, [b for b in cur_blocks if b is not None])

    def output_reference_insertions(cur_chrom, cur_chrom_idx, cur_pos, chrom, pos):
        """Emit blocks for the reference sequence up to the given chrom and pos,
        starting from cur_chrom and cur_pos."""
        # Check if we are in a reference insertion.
        if chrom != chroms[cur_chrom_idx]:
            assert chrom > chroms[cur_chrom_idx]
            # Finish off the chromosome.
            logger.debug('Finishing off chromosome')
            size = chrom_sizes[cur_chrom]
            if cur_pos != size:
                block = Block([BlockLine(genome, cur_chrom, cur_pos, size, '+', 'X' * (size - cur_pos))])
                output.write(str(block))

            # Handle any fully inserted chromosomes.
            cur_chrom_idx += 1
            cur_chrom = chroms[cur_chrom_idx]
            cur_pos = 0
            while cur_chrom != chrom:
                logger.debug('Outputting block for unaligned chrom %s', cur_chrom)
                block = Block([BlockLine(genome, cur_chrom, 0, chrom_sizes[cur_chrom], '+', 'X' * (chrom_sizes[cur_chrom]))])
                output.write(str(block))
                cur_chrom_idx += 1
                cur_chrom = chroms[cur_chrom_idx]
        if pos != cur_pos:
            logger.debug('Outputting block for reference insertion before next block at %s (cur_pos %s)', pos, cur_pos)
            assert pos > cur_pos
            block = Block([BlockLine(genome, cur_chrom, cur_pos, pos, '+', 'X' * (pos - cur_pos))])
            output.write(str(block))

        return cur_chrom, cur_chrom_idx, cur_pos

    chroms = sorted(chrom_sizes.keys())
    cur_chrom_idx = 0
    cur_chrom = chroms[0]
    cur_pos = 0

    # Sentinel indicating that we should fetch the next block from this file.
    # FIXME: this is awful, but shit needs to get done
    need_next = 'need next block'
    # Streams of Blocks from each file, in order.
    block_streams = [merged_dup_stream(f, Block.read_next_from_file) for f in block_files]
    # Our current block from each file, in order (or None if we are at EOF, or
    # need_next if the next block is desired).
    cur_blocks = [next(stream, None) for stream in block_streams]

    while True:
        if all([b is None for b in cur_blocks]):
            # Finished.
            break

        # Find a "next" block after our last position.
        smallest = get_smallest_block(cur_blocks)

        # Handle any insertions before the start of the next block.
        cur_chrom, cur_chrom_idx, cur_pos = output_reference_insertions(cur_chrom, cur_chrom_idx, cur_pos, smallest.first.chrom, smallest.first.start)

        if not any(smallest.overlaps(b) for b in cur_blocks if b is not None and b != smallest):
            # Print out the smallest block if it doesn't overlap anything.
            logger.debug('Found non-overlapping block')
            output.write(str(smallest))
            smallest_idx = cur_blocks.index(smallest)
            cur_blocks[smallest_idx] = need_next
            cur_pos = smallest.first.end
        else:
            # Split up the current blocks at the end of their overlap.
            blocks_to_merge = []
            split_point = find_split_point(cur_blocks)
            assert smallest.first.start < split_point <= smallest.first.end
            for i, block in enumerate(cur_blocks):
                if block is None:
                    # Already done with this file
                    continue
                if smallest.first.overlaps(block.first) and block.first.start < split_point < block.first.end:
                    left_block, right_block = block.split(split_point - block.first.start)
                    blocks_to_merge.append(left_block)
                    # Keep the right side of the block around for the next iteration.
                    cur_blocks[i] = right_block
                elif smallest.first.overlaps(block.first) and split_point == block.first.end:
                    # We shouldn't split this block, because it's already a perfect overlap.
                    blocks_to_merge.append(block)
                    cur_blocks[i] = need_next

            logger.debug('Merging %s overlapping blocks', len(blocks_to_merge))

            # Merge the blocks and emit the result.
            growing_block = blocks_to_merge[0]
            for block in blocks_to_merge[1:]:
                growing_block = growing_block.merge(block)
            cur_pos = growing_block.first.end
            output.write(str(growing_block))

        # Done with this iteration. Advance in whichever of the files we need
        # blocks from. (We may have blocks left over from the splits, so we
        # shouldn't advance in all of them all the time.)
        for i, block_or_sentinel in enumerate(cur_blocks):
            if block_or_sentinel == need_next:
                # Read from the corresponding file
                cur_blocks[i] = next(block_streams[i], None)

    # Handle any insertions after all the alignments are done.
    output_reference_insertions(cur_chrom, cur_chrom_idx, cur_pos, chroms[-1], chrom_sizes[chroms[-1]])

def maximize_block_length(blocks, ref_sequence, output):
    """
    Go through the blocks in the given file and concatenate them together "horizontally" if they are adjacent.

    Also fills in the reference sequence lines with the correct bases from ref_sequence.
    """
    # Our current block. This will grow as we find more blocks to add to it.
    growing_block = None
    while True:
        cur_block = Block.read_next_from_file(blocks)
        if cur_block is None:
            break
        logger.debug('cur_block: %s, growing_block: %s', cur_block, growing_block)
        if growing_block is not None:
            assert cur_block.first.chrom >= growing_block.first.chrom and \
                (cur_block.first.start == growing_block.first.end or cur_block.first.start == 0), \
                "Input blocks must be sorted and span the entire sequence without overlap"
        # Change the reference bases from 'X' to the real sequence
        cur_block.fill_in_first_sequence(ref_sequence[cur_block.first.chrom][cur_block.first.start:cur_block.first.end])
        if growing_block is None:
            growing_block = cur_block
        elif growing_block.compatible_with(cur_block):
            # This block is OK to merge with our previous block, so go ahead and
            # do it.
            growing_block.concatenate(cur_block)
        else:
            # Not fully adjacent to the previous block, so we can't concatenate
            # them. Emit our block and start over.
            output.write(str(growing_block))
            growing_block = cur_block
    output.write(str(growing_block))

def maf_export_job(job, hal, genomes, blocksID, opts):
    """Once blocks are on the reference, export them to MAF."""
    chrom_sizes = {}
    for genome in genomes:
        chrom_sizes[genome] = get_chrom_sizes(hal, genome)

    blocks_path = job.fileStore.readGlobalFile(blocksID)
    output_path = job.fileStore.getLocalTempFile()
    with open(blocks_path) as blocks, open(output_path, 'w') as output:
        output.write('##maf version=1\n')
        while True:
            block = Block.read_next_from_file(blocks)
            if block is None:
                break
            output.write(block.to_maf(chrom_sizes))
    return job.fileStore.writeGlobalFile(output_path)

class Alignment(namedtuple('Alignment', ('my_chrom', 'my_start', 'my_end',
                                         'other_chrom', 'other_start', 'other_end',
                                         'other_strand'))):
    """
    Represents a single pairwise alignment.

    The alignment is always from the positive strand of the reference (i.e. the "my_" genome).
    """
    def split(self, split_point):
        """
        Return two new alignments created by splitting this alignment at an offset relative to its start.
        """
        if self.other_strand == '+':
            return (Alignment(self.my_chrom, self.my_start, self.my_start + split_point,
                              self.other_chrom, self.other_start, self.other_start + split_point,
                              self.other_strand),
                    Alignment(self.my_chrom, self.my_start + split_point, self.my_end,
                              self.other_chrom, self.other_start + split_point, self.other_end,
                              self.other_strand))
        else:
            return (Alignment(self.my_chrom, self.my_start, self.my_start + split_point,
                              self.other_chrom, self.other_end - split_point, self.other_end,
                              self.other_strand),
                    Alignment(self.my_chrom, self.my_start + split_point, self.my_end,
                              self.other_chrom, self.other_start, self.other_end - split_point,
                              self.other_strand))

class AlignmentGroup(object):
    """
    Represents a collection of overlapping pairwise alignments (or just one).

    The overlap is always perfect on the reference, that is, each pairwise
    alignment refers to the exact same reference region.
    """
    def __init__(self, alignments):
        self.alignments = alignments
        self.chrom = alignments[0].my_chrom
        self.start = alignments[0].my_start
        self.end = alignments[0].my_end

    def split(self, split_point):
        """
        Return two new AlignmentGroups created by splitting them at an offset relative to its start.
        """
        splits = [aln.split(split_point) for aln in self.alignments]
        left = [alns[0] for alns in splits]
        right = [alns[1] for alns in splits]
        return AlignmentGroup(left), AlignmentGroup(right)

    def overlaps(self, other):
        """
        Do these two AlignmentGroups share any bases on the reference?
        """
        if self.chrom == other.chrom \
           and ((self.start < other.end <= self.end) \
                or (self.start <= other.start < self.end)):
            return True
        else:
            return False

    def merge(self, other):
        """
        Add the alignments from another AlignmentGroup into this one.

        The alignments must perfectly overlap these alignments on the reference.
        """
        assert self.chrom == other.chrom \
            and self.start == other.start \
            and self.end == other.end
        return AlignmentGroup(self.alignments + other.alignments)

    def __repr__(self):
        return self.__class__.__name__ + "(" + ",".join([repr(a) for a in self.alignments]) + ")"

def read_next_alignment(alignmentFile):
    """
    Get the next AlignmentGroup from alignmentFile.
    """
    line = alignmentFile.readline()
    if len(line) == 0:
        # EOF
        return None
    fields = line.strip().split('\t')
    assert len(fields) == 7, "Misformatted alignment line"
    fields[1], fields[2], fields[4], fields[5] = int(fields[1]), int(fields[2]), int(fields[4]), int(fields[5])
    alignment = Alignment(*fields)
    assert alignment.my_end - alignment.my_start == alignment.other_end - alignment.other_start
    return AlignmentGroup([alignment])

@attrs(hash=True, slots=True)
class BlockLine(object):
    """
    Represents a single sequence line within a block.

    This class avoids copying the input sequence on splits, to avoid quadratic
    behavior when a large block is repeatedly split.
    """
    genome = attrib()
    chrom = attrib()
    start = attrib()
    end = attrib()
    strand = attrib()
    # Aligned sequence for this block entry.
    aligned_seq = attrib(cmp=False)
    # Our window within the above aligned sequence. This allows us to keep
    # references to the same string within many different block entries, rather
    # than copying giant strings every time we split or reverse. Quadratic
    # behavior would really hurt on chromosomes that are hundreds of megabases
    # long.
    align_start = attrib(default=0, cmp=False)
    align_end = attrib(cmp=False)

    @align_end.default
    def _align_end(self):
        return len(self.aligned_seq)

    # Whether we need to be reverse complemented eventually. This is computed
    # lazily to avoid quadratic behavior.
    needs_rev_comp = attrib(default=False)

    def seq_pos_to_align_pos(self, seq_pos):
        """Given a position in our chromosome, get the corresponding offset within our aligned sequence."""
        align_pos = self.align_start
        if self.strand == '+':
            my_seq_pos = self.start
            step = +1
        else:
            my_seq_pos = self.end - 1
            step = -1
        while my_seq_pos != seq_pos:
            if self.aligned_seq[align_pos] != '-':
                my_seq_pos += step
            align_pos += 1

        return align_pos - self.align_start

    def reverse(self):
        """Return a new BlockLine referring to the opposite strand."""
        return BlockLine(
            self.genome,
            self.chrom,
            self.start,
            self.end,
            '+' if self.strand == '-' else '-',
            self.aligned_seq,
            align_start=self.align_start,
            align_end=self.align_end,
            needs_rev_comp=True
        )

    def align_pos_to_seq_pos(self, tgt_align_pos):
        """Given a offset within our aligned sequence, get the corresponding position in our chromosome."""
        if self.strand == '+':
            seq_pos = self.start
        else:
            seq_pos = self.end
        align_pos = 0
        while align_pos < tgt_align_pos:
            if self.aligned_seq[self.align_start + align_pos] != '-':
                if self.strand == '+':
                    seq_pos += 1
                else:
                    seq_pos -= 1
            align_pos += 1
        return seq_pos

    def split(self, split_pos):
        """
        Return two new BlockLines referring to the alignment before and after split_pos.

        split_pos refers to an offset *IN ALIGNED CHARACTERS*, not in sequence
        positions, from the start position of this block.
        """
        # Find position (relative to sequence, not alignment) in this line at split point
        my_pos = self.align_pos_to_seq_pos(split_pos)
        # Create 2 new block lines
        if self.strand == '+':
            first_block_line = BlockLine(self.genome, self.chrom, self.start, my_pos, '+', self.aligned_seq, align_start=self.align_start, align_end=self.align_start + split_pos)
            second_block_line = BlockLine(self.genome, self.chrom, my_pos, self.end, '+', self.aligned_seq, align_start=self.align_start + split_pos, align_end=self.align_end)
        else:
            first_block_line = BlockLine(self.genome, self.chrom, my_pos, self.end, '-', self.aligned_seq, align_start=self.align_start, align_end=self.align_start + split_pos)
            second_block_line = BlockLine(self.genome, self.chrom, self.start, my_pos, '-', self.aligned_seq, align_start=self.align_start + split_pos, align_end=self.align_end)
        return first_block_line, second_block_line

    def overlaps(self, other):
        """Does this BlockLine share any bases with the other?"""
        if self.genome == other.genome and self.chrom == other.chrom \
           and ((self.start < other.end <= self.end) \
                or (self.start <= other.start < self.end)):
            return True
        else:
            return False

    def precedes(self, other):
        """Does this BlockLine *immediately* precede the other?"""
        if self.genome == other.genome \
           and self.chrom == other.chrom \
           and self.strand == other.strand:
            if self.strand == '+':
                return other.start == self.end
            else:
                return other.end == self.start
        else:
            return False

    @property
    def seq(self):
        """Get the actual aligned characters this BlockLine refers to."""
        seq = self.aligned_seq[self.align_start:self.align_end]
        if self.needs_rev_comp:
            return reverseComplement(seq)
        else:
            return seq

    @seq.setter
    def seq(self, seq):
        self.aligned_seq = seq
        self.align_start = 0
        self.align_end = len(seq)
        self.needs_rev_comp = False

@attrs
class Block(object):
    block_lines = attrib()

    @classmethod
    def read_next_from_file(cls, f):
        """Read the next block stanza from f into memory."""
        block_lines = []
        while True:
            line = f.readline()
            if len(line) == 0:
                # EOF
                return None
            if len(line.strip()) == 0:
                # End of block stanza
                break
            fields = line.strip().split('\t')
            fields[2], fields[3] = int(fields[2]), int(fields[3])
            block_lines.append(BlockLine(genome=fields[0],
                                         chrom=fields[1],
                                         start=fields[2],
                                         end=fields[3],
                                         strand=fields[4],
                                         aligned_seq=fields[5]))

        # Sanity check
        seq_len = len(block_lines[0].seq)
        assert all(len(i.seq) == seq_len for i in block_lines)

        return cls(block_lines)

    @property
    def first(self):
        return self.block_lines[0]

    @property
    def chrom(self):
        return self.first.chrom

    @property
    def start(self):
        return self.first.start

    @property
    def end(self):
        return self.first.end

    def split(self, pos):
        """
        Split at a given index *in sequence positions* relative to the first block entry.

        Creates two new blocks, one of which is everything left of the split
        point, non-inclusive, and the other of which is everything right of the
        split point, inclusive.
        """
        assert self.first.strand == '+'
        assert 0 <= pos < self.first.end - self.first.start, "Split point must be within block"
        align_pos = self.first.seq_pos_to_align_pos(self.first.start + pos)
        first_block_lines = []
        second_block_lines = []
        for block_line in self.block_lines:
            first_line, second_line = block_line.split(align_pos)
            first_block_lines.append(first_line)
            second_block_lines.append(second_line)

        # Remove any completely empty lines (signified by start == end). These
        # appear when all the non-gap characters of a line are right or left of
        # the split point.
        first_block_lines = [bl for bl in first_block_lines if bl.start != bl.end]
        second_block_lines = [bl for bl in second_block_lines if bl.start != bl.end]

        return Block(first_block_lines), Block(second_block_lines)

    def prepend_new_sequence(self, ref_start, ref_end, new_genome, new_chrom, new_start, new_end, strand):
        """
        Add a new reference sequence line to this block.

        The new reference sequence will not have any real bases, just 'X' and
        gap characters. The new sequence is aligned to the old reference
        sequence starting at ref_start and ending at ref_end (in the orientation
        given by the strand parameter).
        """
        assert ref_end - ref_start == new_end - new_start
        assert 0 <= ref_start < ref_end
        assert 0 <= new_start < new_end
        pre_gap = self.first.seq_pos_to_align_pos(ref_start)
        ref_end_align_pos = self.first.seq_pos_to_align_pos(ref_end)
        post_gap = len(self.first.seq) - ref_end_align_pos
        block_lines = [bl for bl in self.block_lines]
        aligned_portion = re.sub(r'[^-]', 'X', self.first.seq[pre_gap:ref_end_align_pos])
        assert len(aligned_portion) + pre_gap + post_gap == len(self.first.seq), \
            "Prepended sequence must match the current length of the block!"
        if strand == '-':
            # Flip everything. The first sequence should always be in the + orientation
            block_lines = [bl.reverse() for bl in block_lines]
            new_seq = '-' * post_gap + "".join(list(reversed(aligned_portion))) + '-' * pre_gap
        else:
            new_seq = '-' * pre_gap + aligned_portion + '-' * post_gap
        block_lines.insert(0, BlockLine(new_genome, new_chrom, new_start, new_end, '+', new_seq))
        new_block = Block(block_lines)
        return new_block

    def __str__(self):
        s = ""
        for block_line in self.block_lines:
            assert len(block_line.seq) == len(self.first.seq)
            s += "{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{seq}\n".format(
                genome=block_line.genome, chrom=block_line.chrom, start=block_line.start,
                end=block_line.end, strand=block_line.strand, seq=block_line.seq)
        s += "\n"
        return s

    def fill_in_first_sequence(self, first_seq):
        """Set the non-gap bases of the reference line to be equal to the sequence in first_seq."""
        assert len(first_seq) == self.first.end - self.first.start
        assert self.first.strand == '+'
        aligned_seq = self.first.seq
        # Convert to list and assign the correct characters, skipping gaps
        new_aligned_seq = list(aligned_seq)
        seq_pos = 0
        for i, align_char in enumerate(new_aligned_seq):
            if align_char != '-':
                new_aligned_seq[i] = first_seq[seq_pos]
                seq_pos += 1
        # Convert back to a string
        new_aligned_seq = "".join(new_aligned_seq)
        self.block_lines[0] = BlockLine(self.first.genome, self.first.chrom, self.first.start, self.first.end, self.first.strand, new_aligned_seq)

    def merge(self, other):
        """
        Return a new block adding in the alignments from another block.

        The other block must overlap with or abut this block (in the
        reference). Its alignments are stacked "vertically" after this block's
        alignments in the new block. The overlap doesn't need to be perfect; gap
        characters to maintain alignment relationships where necessary.
        """
        # Ensure there's some overlap or adjacency
        assert (self.first.start <= other.first.end <= self.first.end) \
            or (self.first.start <= other.first.start <= self.first.end), \
            "Blocks must overlap or abut on the reference"
        assert self.first.genome == other.first.genome, \
            "Blocks must be referenced on the same genome"
        assert self.first.chrom == other.first.chrom, \
            "Blocks must be on the same reference chromosome"
        assert self.first.strand == '+' and other.first.strand == '+', \
            "Blocks must have their first line in the + orientation"

        # Ensure there's no overlap between non-reference lines. We're
        # blindly stacking these lines on top of each other, so we'll
        # create an invalid block if there is indeed overlap.
        for my_bl in self.block_lines[1:]:
            for other_bl in other.block_lines[1:]:
                assert not my_bl.overlaps(other_bl)

        new_start = min(self.first.start, other.first.start)
        new_end = max(self.first.end, other.first.end)
        # To avoid ugly quadratic behavior we represent the new sequences as
        # lists and join them into strings only at the end.
        block_seqs = [[] for _ in range(len(self.block_lines) + len(other.block_lines) - 1)]
        my_seq_pos = self.first.start
        other_seq_pos = other.first.start
        my_align_pos = 0
        other_align_pos = 0
        my_end_pos = len(self.first.seq)
        other_end_pos = len(other.first.seq)

        # Region before overlap.
        overlap_start = max(self.first.start, other.first.start)
        # TODO: think about how to unify this logic (as well as that for the
        # post-overlap region below); it's difficult to abstract this, but even
        # more difficult to tell if there are subtle bugs. We could choose to
        # handle only the case where self.first.start < other.first.start,
        # swapping the merger and mergee if not, but that would only get rid of
        # one branch here.
        if my_seq_pos < overlap_start:
            pre_overlap_align_len = self.first.seq_pos_to_align_pos(overlap_start)
            # First block's sequences: taken from the first block's alignment
            for i in range(len(self.block_lines)):
                block_seqs[i].extend(self.block_lines[i].seq[:pre_overlap_align_len])
            # Second block's sequences: all gaps
            for i in range(len(self.block_lines), len(block_seqs)):
                block_seqs[i].extend('-' * pre_overlap_align_len)
            my_align_pos += pre_overlap_align_len
            my_seq_pos = overlap_start
        elif other_seq_pos < overlap_start:
            pre_overlap_align_len = other.first.seq_pos_to_align_pos(overlap_start)
            # Reference sequence: taken from the second block
            block_seqs[0].extend(other.block_lines[0].seq[:pre_overlap_align_len])
            # First block's sequences: all gaps
            for i in range(1, len(self.block_lines)):
                block_seqs[i].extend('-' * pre_overlap_align_len)
            # Second block's sequences: taken from the second block's alignment
            for i in range(len(self.block_lines), len(block_seqs)):
                block_seqs[i].extend(other.block_lines[i - len(self.block_lines) + 1].seq[:pre_overlap_align_len])
            other_align_pos += pre_overlap_align_len
            other_seq_pos = overlap_start

        # Region of overlap.

        def update_columns_from_reference_deletions(my_align_pos, other_align_pos, block_seqs):
            # Handle reference deletions in first block
            while my_align_pos < my_end_pos and self.block_lines[0].seq[my_align_pos] == '-':
                column = [bl.seq[my_align_pos] for bl in self.block_lines]
                for i, char in enumerate(column + ['-'] * (len(other.block_lines) - 1)):
                    block_seqs[i].append(char)
                my_align_pos += 1

            # Handle reference deletions within second block
            while other_align_pos < other_end_pos and other.block_lines[0].seq[other_align_pos] == '-':
                column = [bl.seq[other_align_pos] for bl in other.block_lines]
                for i, char in enumerate(['-'] * (len(self.block_lines) - 1) + column):
                    block_seqs[i].append(char)
                other_align_pos += 1
            return my_align_pos, other_align_pos

        my_align_pos, other_align_pos = update_columns_from_reference_deletions(my_align_pos, other_align_pos, block_seqs)

        overlap_end = min(self.first.end, other.first.end)
        while my_seq_pos < overlap_end:
            assert other_seq_pos < overlap_end
            # Add the current column from both blocks, then all insertions
            # relative to the reference in the first block, then all insertions
            # relative to the reference in the second block.
            my_cur_column = [bl.seq[my_align_pos] for bl in self.block_lines]
            assert my_cur_column[0] != '-'
            other_cur_column = [bl.seq[other_align_pos] for bl in other.block_lines]
            assert other_cur_column[0] != '-'
            assert my_cur_column[0] == other_cur_column[0], "Mismatch in reference base when merging blocks"
            for i, char in enumerate(my_cur_column + other_cur_column[1:]):
                block_seqs[i].append(char)

            my_align_pos += 1
            other_align_pos += 1
            my_align_pos, other_align_pos = update_columns_from_reference_deletions(my_align_pos, other_align_pos, block_seqs)
            my_seq_pos += 1
            other_seq_pos += 1

        # Region after overlap.
        if my_align_pos < my_end_pos:
            # First block's sequences: taken from the first block's alignment
            for i in range(len(self.block_lines)):
                block_seqs[i].extend(self.block_lines[i].seq[my_align_pos:])
            # Second block's sequences: all gaps
            for i in range(len(self.block_lines), len(block_seqs)):
                block_seqs[i].extend('-' * (my_end_pos - my_align_pos))
        elif other_align_pos < other_end_pos:
            # Reference sequence: taken from the second block
            block_seqs[0].extend(other.block_lines[0].seq[other_align_pos:])
            # First block's sequences: all gaps
            for i in range(1, len(self.block_lines)):
                block_seqs[i].extend('-' * (other_end_pos - other_align_pos))
            # Second block's sequences: taken from the second block's alignment
            for i in range(len(self.block_lines), len(block_seqs)):
                block_seqs[i].extend(other.block_lines[i - len(self.block_lines) + 1].seq[other_align_pos:])

        # Build up and return the new block.
        block_seqs = ["".join(seq) for seq in block_seqs]

        ref_block_line = BlockLine(self.first.genome, self.first.chrom, new_start, new_end, '+', block_seqs[0])
        new_block_lines = [ref_block_line] + \
                          [BlockLine(bl.genome, bl.chrom, bl.start, bl.end, bl.strand, block_seqs[i + 1]) for i, bl in enumerate(self.block_lines[1:])] + \
                          [BlockLine(bl.genome, bl.chrom, bl.start, bl.end, bl.strand, block_seqs[i + len(self.block_lines)]) for i, bl in enumerate(other.block_lines[1:])]
        return Block(new_block_lines)

    def compatible_with(self, other):
        """
        Return whether this block can be concatened with the other.
        """
        # Our reference entry absolutely has to be syntenic with the other.
        if self.overlaps(other) or not self.first.precedes(other.first):
            return False

        # Gather up all our block lines, indexed by their genome.
        genome_to_block_lines = defaultdict(list)
        for block_line in self.block_lines:
            genome_to_block_lines[block_line.genome].append(block_line)

        # Now go through all the other block's lines, verifying that they can
        # follow the corresponding line (if any) in our block.
        for block_line in other.block_lines:
            potential_partner_lines = genome_to_block_lines[block_line.genome]
            partner = None
            for potential_partner in potential_partner_lines:
                if block_line.overlaps(potential_partner):
                    # Can't have two copies of the same base in a block.
                    return False
                if potential_partner.precedes(block_line):
                    if partner is None or partner.precedes(potential_partner):
                        partner = potential_partner
            # Check for a rearrangement relative to a previous block. Though we
            # could stack these and still create a valid block, it's probably
            # not a great idea.
            if len(potential_partner_lines) != 0 and partner is None:
                return False
        return True

    def overlaps(self, other):
        """
        Do these blocks' reference lines overlap?
        """
        if self.first.overlaps(other.first):
            return True
        else:
            return False

    def find_partners(self, other):
        corresponding_lines = {}
        # Gather up all our block lines, indexed by their genome.
        genome_to_block_lines = defaultdict(list)
        for block_line in self.block_lines:
            genome_to_block_lines[block_line.genome].append(block_line)

        for block_line in other.block_lines:
            potential_partner_lines = genome_to_block_lines[block_line.genome]
            partner = None
            for potential_partner in potential_partner_lines:
                if potential_partner.precedes(block_line):
                    assert partner is None
                    partner = potential_partner
            if partner is not None:
                corresponding_lines[block_line] = partner

        return corresponding_lines

    def concatenate(self, other):
        """
        Combine the other block's lines with ours, extending our lines to include theirs.

        The other block must abut this one in any matching lines. Any lines that
        don't have a match between the two blocks will be added as well, but
        gapped out in the positions where they are missing.
        """
        assert self.first.end == other.first.start
        corresponding_lines = self.find_partners(other)
        lines_to_add = []
        my_block_length = len(self.first.seq)
        other_block_length = len(other.first.seq)
        for other_block_line in other.block_lines:
            if other_block_line in corresponding_lines:
                # We can extend our existing block line with this one.
                self_block_line = corresponding_lines[other_block_line]
                if self_block_line.strand == '+':
                    self_block_line.end = other_block_line.end
                else:
                    self_block_line.start = other_block_line.start
                self_block_line.seq += other_block_line.seq
            else:
                # No partner line; we just stack it on vertically, taking into
                # account the extra gaps at the beginning.
                other_block_line.seq = '-' * (my_block_length) + other_block_line.seq
                lines_to_add.append(other_block_line)

        # Go through and add extra gaps to any of our lines that weren't
        # partnered with any from the other block.
        self_lines_with_partners = set(corresponding_lines.values())
        for self_block_line in self.block_lines:
            if self_block_line not in self_lines_with_partners:
                self_block_line.seq = self_block_line.seq + '-' * (other_block_length)

        self.start_pos = 0
        self.end_pos = len(self_block_line.seq)

        for line in lines_to_add:
            self.block_lines.append(line)

    def to_maf(self, chrom_sizes):
        """Convert this block into a MAF block, getting the chrom length from the chrom_sizes dict."""
        ret = "a\n"
        for block_line in self.block_lines:
            chr_size = chrom_sizes[block_line.genome][block_line.chrom]
            if block_line.strand == '+':
                start = block_line.start
            else:
                start = chr_size - block_line.end
            ret += "s {genome}.{chr} {start} {size} {strand} {chr_size} {seq}\n".format(
                genome=block_line.genome,
                chr=block_line.chrom,
                start=start,
                size=block_line.end - block_line.start,
                strand=block_line.strand,
                chr_size=chr_size,
                seq=block_line.seq)
        ret += "\n"
        return ret

def lift_region(alignment, ref_start, ref_end):
    """Lift a given region on the reference using the provided Alignment."""
    assert alignment.my_start <= ref_start < alignment.my_end
    assert alignment.my_start < ref_end <= alignment.my_end
    start_offset = ref_start - alignment.my_start
    end_offset = alignment.my_end - ref_end
    if alignment.other_strand == '+':
        other_start = alignment.other_start + start_offset
        other_end = alignment.other_end - end_offset
    else:
        other_start = alignment.other_start + end_offset
        other_end = alignment.other_end - start_offset
    return other_start, other_end

def paired_iter(iterable):
    """Given an iterable of items, return an iterable that gives (item1, item2), (item2, item3)...

    The last pair will be (final_item, None)."""
    first_iter, second_iter = itertools.tee(iterable, 2)
    next(second_iter, None)
    return itertools.zip_longest(first_iter, second_iter)

def lift_blocks(alignments, blocks, other_genome, output):
    """
    Lift the given blocks file to a new reference using the given pairwise alignments file.

    The pairwise alignments must point from the old reference to the new.
    """
    aln_stream = paired_iter(merged_dup_stream(alignments, read_next_alignment))
    aln, next_aln = next(aln_stream, (None, None))
    block = Block.read_next_from_file(blocks)
    while block is not None and aln is not None:
        assert block.first.strand == '+'

        # Make sure we are at a point where we have some overlap between our
        # alignment and our block, advancing through the files and discarding
        # alignments/blocks if necessary.
        need_restart = False
        while block.first.chrom < aln.chrom or (aln.chrom == block.first.chrom and block.first.end <= aln.start):
            # Block doesn't fit within alignment, fast-forward on block side
            block = Block.read_next_from_file(blocks)
            need_restart = True
            if block is None:
                return
        while aln.chrom < block.first.chrom or (aln.chrom == block.first.chrom and aln.end <= block.first.start):
            # Need to fast-forward on alignment side
            aln, next_aln = next(aln_stream, (None, None))
            need_restart = True
            if aln is None:
                return

        if need_restart:
            # We need to go back to the beginning, because after advancing, our
            # alignment may now be past our block, or vice versa.
            continue

        # OK, now we have an alignment which overlaps a block. We can do the actual lifting work below.

        assert aln.chrom == block.chrom
        assert aln.start < block.end
        assert block.first.start < aln.end

        ref_start = max(aln.start, block.start)
        ref_end = min(aln.end, block.end)
        if block.end > aln.end and (next_aln is not None and next_aln.chrom == aln.chrom and next_aln.start < block.end):
            # Block goes over the edge of alignment, so we need to split it
            assert block.end > aln.end
            # We should split based on the *next* alignment. If we
            # have a gap between alignments, i.e. an insertion
            # relative to our target, then we should always split
            # *after* the insertion, to ensure we lift as much
            # syntenic sequence as possible.
            left_block, right_block = block.split(next_aln.start - block.first.start)
            block = left_block
            next_block = right_block
        else:
            # No need to split this block. Either it lies entirely
            # within this alignment, or it goes past the end of the
            # current alignment, but there's no following alignment
            # that we would need to lift a different part of the block
            # to.
            next_block = Block.read_next_from_file(blocks)
        for sub_aln in aln.alignments:
            # For each (possibly) duplicated alignment in this region,
            # output a new block with the reference sequence added.
            other_start, other_end = lift_region(sub_aln, ref_start, ref_end)
            new_block = block.prepend_new_sequence(ref_start, ref_end, other_genome, sub_aln.other_chrom, other_start, other_end, sub_aln.other_strand)
            logger.debug('outputting block %s', new_block)
            output.write(str(new_block))
        block = next_block

def get_sequence(hal_path, genome):
    """
    Get a dict from {seq_name: sequence} representing the nucleotide sequence for a genome.
    """
    with get_ascii_tmp_file() as tmp:
        call(['hal2fasta', hal_path, genome, '--outFaPath', tmp.name])
        return dict(fastaRead(tmp.name))

def get_chrom_sizes(hal_path, genome):
    """Get a dict of {chrom_name: length} for a genome from a hal."""
    bed = call(['halStats', '--bedSequences', genome, hal_path])
    lengths = {}
    for line in bed.split("\n"):
        if len(line.strip()) == 0:
            # Empty line
            continue
        fields = line.split()
        lengths[fields[0]] = int(fields[2])
    return lengths

def get_leaf_blocks(hal_path, genome, output_path):
    """Get the initial blocks file for a leaf genome."""
    seqs = get_sequence(hal_path, genome)
    lengths = get_chrom_sizes(hal_path, genome)
    with open(output_path, 'w') as out:
        for chrom in sorted(lengths.keys()):
            length = lengths[chrom]
            # Write the block stanza. In this case, the block is just the whole chromosome.
            seq = seqs[chrom]
            assert len(seq) == length, "chrom.sizes and fasta must agree on sequence length"
            out.write("{genome}\t{chrom}\t{start}\t{end}\t{strand}\t".format(genome=genome, chrom=chrom, start=0, end=length, strand='+'))
            out.write(seq)
            out.write("\n")
            out.write("\n")

def sort_blocks(blocks, output):
    """Sort a blocks file relative to the first line.

    This uses UNIX sort, because this can potentially need to sort hundreds of
    millions to billions of blocks, which can take >30GB of memory if done in
    Python naively. sort will intelligently use memory/disk resources to avoid
    memory exhaustion.
    """
    # Gather up the start positions within all the blocks, and their position within the file.
    with get_ascii_tmp_file(delete=False) as file_to_sort:
        while True:
            file_pos = blocks.tell()
            block = Block.read_next_from_file(blocks)
            if block is None:
                break
            file_to_sort.write("%s\t%s\t%s\t%s\n" % (block.first.chrom, block.first.start,
                                                     block.first.end, file_pos))
    # Do the actual sorting
    call(['sort', '-k1,1', '-k2,2n', '-k3,3n', file_to_sort.name, '-o', file_to_sort.name])

    with open(file_to_sort.name, encoding='ascii') as sorted:
        for line in sorted:
            fields = line.strip().split()
            file_pos = int(fields[-1])
            blocks.seek(file_pos)
            block = Block.read_next_from_file(blocks)
            output.write(str(block))
    os.remove(file_to_sort.name)

def get_alignment_to_parent(hal_path, child_genome, output_path):
    """Get a file representing maximal gapless alignment blocks between this child and its parent (referenced on the child)."""
    with get_ascii_tmp_file() as tmp:
        call(['halAlignedExtract', '--viewParentCoords', '--alignedFile', tmp.name, hal_path, child_genome])
        with open(output_path, 'w') as output:
            maximize_gapless_alignment_length(open(tmp.name), output)
        call(['sort', output_path, '-t', '\t', '-k1,1', '-k2,2n', '-o', output_path])

def get_alignment_to_child(hal_path, child_genome, output_path):
    """Get a file representing maximal gapless alignment blocks between this child and its parent (referenced on the parent)."""
    with get_ascii_tmp_file() as original_alignment, get_ascii_tmp_file() as flipped_alignment:
        get_alignment_to_parent(hal_path, child_genome, original_alignment.name)
        flip_alignment(open(original_alignment.name), open(flipped_alignment.name, 'w'))
        call(['sort', flipped_alignment.name, '-k1,1', '-k2,2n', '-o', output_path])

def flip_alignment(input, output):
    """
    Flip alignment from halAlignedExtract so that it is parent->child rather than child->parent.
    """
    for line in input:
        line = line.strip()
        if len(line) == 0:
            # Blank line
            continue
        fields = line.split('\t')
        child_chr, child_start, child_end, \
            parent_chr, parent_start, parent_end, parent_strand = fields
        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (parent_chr, parent_start, parent_end, child_chr, child_start, child_end, parent_strand))

def maximize_gapless_alignment_length(input, output):
    """Merge fully-adjacent blocks from halAlignedExtract."""
    cur_block = None
    for line in input:
        line = line.strip()
        if len(line) == 0:
            # Blank line
            continue
        fields = line.split('\t')
        new_child_chr, new_child_start, new_child_end, \
            new_parent_chr, new_parent_start, new_parent_end, new_parent_strand = fields
        if cur_block is None:
            cur_block = fields
            continue
        cur_child_chr, cur_child_start, cur_child_end, \
            cur_parent_chr, cur_parent_start, cur_parent_end, cur_parent_strand = cur_block
        # Can we merge this new block with our current block?
        if cur_child_chr == new_child_chr \
           and cur_child_end == new_child_start \
           and cur_parent_chr == new_parent_chr \
           and cur_parent_end == new_parent_start \
           and cur_parent_strand == new_parent_strand:
            # Mergeable.
            cur_block = (cur_child_chr, cur_child_start, new_child_end, cur_parent_chr, cur_parent_start, new_parent_end, cur_parent_strand)
        else:
            # Not mergeable: print the block we finished extending, and prepare this block to be the new block.
            output.write("\t".join(cur_block) + "\n")
            cur_block = fields
    output.write("\t".join(cur_block) + "\n")

def makeURL(path_or_url):
    if urlparse(path_or_url).scheme == '':
        return "file://" + os.path.abspath(path_or_url)
    else:
        return path_or_url

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('hal')
    parser.add_argument('refGenome')
    parser.add_argument('outputMaf')
    parser.add_argument('--intermediateResultsUrl')
    parser.add_argument('--targetGenomes')
    Job.Runner.addToilOptions(parser)
    opts = parser.parse_args()
    if opts.targetGenomes is not None:
        opts.targetGenomes = opts.targetGenomes.split(",")
    with Toil(opts) as toil:
        if opts.restart:
            mafID = toil.restart()
        else:
            mafID = toil.start(Job.wrapJobFn(start_job, opts.hal, opts.refGenome, opts))
        toil.exportFile(mafID, makeURL(opts.outputMaf))

# Tests

def cleanup_tree(tree):
    """Make the tree's newick string predictable by sorting children and intifying branch lengths."""
    nodes = list(tree.walk())
    for node in nodes:
        if node.length:
            node.length = int(node.length)
        node.descendants.sort(key=lambda n: n.name)

def test_reroot_tree():
    tree = newick.loads("((C:1,D:1)B:2,(F:3,G:3,(I:3,(K:4)J:2)H:5)E:6)A:9;")[0]
    rerooted = reroot_tree(tree, "J")
    # Ensure the original tree isn't changed.
    assert newick.dumps(tree) == "((C:1,D:1)B:2,(F:3,G:3,(I:3,(K:4)J:2)H:5)E:6)A:9;"
    cleanup_tree(rerooted)
    # Make sure the reroot occurred properly.
    assert newick.dumps(rerooted) == "(((((C:1,D:1)B:2)A:6,F:3,G:3)E:5,I:3)H:2,K:4)J;"

def test_find_node_by_name():
    tree = newick.loads("(A, B, (C, D)E)F;")[0]

    root = find_node_by_name(tree, 'F')
    assert root.name == 'F'
    assert root.ancestor is None

    c = find_node_by_name(tree, 'C')
    assert c.name == 'C'
    assert c.ancestor.name == 'E'

    e = find_node_by_name(tree, 'E')
    assert e.name == 'E'
    assert c in e.descendants
    assert e.ancestor == root

def test_maximize_gapless_alignment_length():
    input = StringIO(dedent("""
    simHuman.chr6	601632	601633	Anc1refChr1	332920	332921	+
    simHuman.chr6	601633	601639	Anc1refChr1	332921	332927	+
    simHuman.chr6	601639	601664	Anc1refChr1	332927	332952	-
    simHuman.chr6	601664	601670	Anc1refChr1	332952	332958	-
    simHuman.chr6	601670	601674	Anc1refChr1	332958	332962	+
    simHuman.chr6	601674	601684	Anc1refChr1	332963	332973	+
    simHuman.chr6	601685	601696	Anc1refChr1	332973	332984	+
    simHuman.chr6	601696	601703	Anc1refChr1	332984	332991	+
    simHuman.chr6	601703	601714	Anc1refChr1	333002	333013	+
    simHuman.chr6	601715	601827	Anc1refChr1	333013	333125	+
    """))
    output = StringIO()

    maximize_gapless_alignment_length(input, output)
    assert output.getvalue() == dedent("""
    simHuman.chr6	601632	601639	Anc1refChr1	332920	332927	+
    simHuman.chr6	601639	601670	Anc1refChr1	332927	332958	-
    simHuman.chr6	601670	601674	Anc1refChr1	332958	332962	+
    simHuman.chr6	601674	601684	Anc1refChr1	332963	332973	+
    simHuman.chr6	601685	601703	Anc1refChr1	332973	332991	+
    simHuman.chr6	601703	601714	Anc1refChr1	333002	333013	+
    simHuman.chr6	601715	601827	Anc1refChr1	333013	333125	+
    """).lstrip()

def test_block_line_split():
    # positive strand
    # start:   v
    # split:   |                v
    seq = 'CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG'
    block_line = BlockLine('Human', 'chr6', 50, 100, '+', seq, align_start=4)
    first, second = block_line.split(17)
    assert first.genome == second.genome == 'Human'
    assert first.chrom == second.chrom == 'chr6'
    assert first.strand == second.strand == '+'
    # Ensure the string is passed around rather than being copied.
    assert first.aligned_seq is second.aligned_seq is seq
    assert first.start == 50
    assert first.end == 61
    assert second.start == 61
    assert second.end == 100

    # Test that the align_start and align_end combine to create the correct seq
    # subsets
    assert first.seq == '--GA---TCAACCAAG-'
    assert second.seq == '-AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG'

    # negative strand
    block_line = BlockLine('Human', 'chr6', 50, 100, '-', seq, align_start=4)
    first, second = block_line.split(17)
    assert first.genome == second.genome == 'Human'
    assert first.chrom == second.chrom == 'chr6'
    assert first.strand == second.strand == '-'
    # Ensure the string is passed around rather than being copied.
    assert first.aligned_seq is second.aligned_seq is seq
    assert first.start == 89
    assert first.end == 100
    assert second.start == 50
    assert second.end == 89

    # Test that the align_start and align_end combine to create the correct seq
    # subsets
    assert first.seq == '--GA---TCAACCAAG-'
    assert second.seq == '-AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG'

def test_block_split():
    input = StringIO(dedent("""
    Human	chr6	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    """).lstrip())

    block = Block.read_next_from_file(input)
    left_block, right_block = block.split(9)
    assert str(left_block) == dedent("""
    Human	chr6	50	59	+	CGG---GA---TCAA
    Chimp	chr6	7045	7059	-	CGGTTGGACCCTCA-

    """).lstrip()
    assert str(right_block) == dedent("""
    Human	chr6	59	100	+	CCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7045	-	--AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    """).lstrip()

    # Now try splitting the right block, which has a start_pos != 0
    block3, block4 = right_block.split(31)
    assert str(block3) == dedent("""
    Human	chr6	59	90	+	CCAAG--AAGAATTGCTATAcctccTAATCCTA
    Chimp	chr6	7017	7045	-	--AAGttAAGAATT---ATAcctccTAATCCTA

    """).lstrip()
    assert str(block4) == dedent("""
    Human	chr6	90	100	+	---TAAAGT----AGCG
    Chimp	chr6	7000	7017	-	cccTAAAGTgggtAGCG

    """).lstrip()

def test_block_prepend_new_sequence():
    # Reverse strand
    input = StringIO(dedent("""
    Human	chr6	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    """).lstrip())

    block = Block.read_next_from_file(input)
    block = block.prepend_new_sequence(56, 61, 'Mouse', 'chrX', 50, 55, '-')
    assert str(block) == dedent("""
    Mouse	chrX	50	55	+	------------------------------------------------XXXXX------------
    Human	chr6	50	100	-	CGCT----ACTTTA---TAGGATTAggaggTATAGCAATTCTT--CTTGGTTGA---TC---CCG
    Chimp	chr6	7000	7059	+	CGCTacccACTTTAgggTAGGATTAggaggTAT---AATTCTTaaCTT---TGAGGGTCCAACCG

    """).lstrip()

    # Positive strand, intersecting with gaps in the reference sequence
    input_2 = StringIO(dedent("""
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    """).lstrip())
    block_2 = Block.read_next_from_file(input_2)
    block_2 = block_2.prepend_new_sequence(0, 8, 'Mouse', 'chrX', 50, 58, '+')
    assert str(block_2) == dedent("""
    Mouse	chrX	50	58	+	X------XX-XXXXX
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    """).lstrip()

def test_block_line_reverse():
    seq = 'GA---TCAACCAAG--AAG'
    block_line = BlockLine('Human', 'chr6', 50, 100, '+', seq, align_start=4, align_end=15)
    reversed = block_line.reverse()
    assert reversed.strand == '-'
    assert reversed.seq == '-CTTGGTTGA-'

def test_prune_tree():
    tree = newick.loads('((A,B)C, (D, E)F, (G)H)I;')[0]

    prune_tree(tree, ['A', 'G'])
    assert newick.dumps(tree) == '((A)C,(G)H)I;'

def test_sort_blocks():
    input = StringIO(dedent("""
    Human	chrX	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    Human	chr6	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    Human	chr6	50	99	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGC
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGC

    Human	chr6	0	50	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    """).lstrip())

    output = StringIO()
    sort_blocks(input, output)

    assert output.getvalue() == dedent("""
    Human	chr6	0	50	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    Human	chr6	50	99	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGC
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGC

    Human	chr6	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    Human	chrX	50	100	+	CGG---GA---TCAACCAAG--AAGAATTGCTATAcctccTAATCCTA---TAAAGT----AGCG
    Chimp	chr6	7000	7059	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcctccTAATCCTAcccTAAAGTgggtAGCG

    """).lstrip()

def test_block_fill_in_first_sequence():
    input = StringIO(dedent("""
    Mouse	chrX	50	55	+	------------------------------------------------X-XXXX------------
    Human	chr6	50	100	-	CGCT----ACTTTA---TAGGATTAggaggTATAGCAATTCTT--C-TTGGTTGA---TC---CCG
    Chimp	chr6	7000	7059	+	CGCTacccACTTTAgggTAGGATTAggaggTAT---AATTCTTaa-CTT---TGAGGGTCCAACCG

    """).lstrip())
    block = Block.read_next_from_file(input)
    block.fill_in_first_sequence("AGTCA")
    assert str(block) == dedent("""
    Mouse	chrX	50	55	+	------------------------------------------------A-GTCA------------
    Human	chr6	50	100	-	CGCT----ACTTTA---TAGGATTAggaggTATAGCAATTCTT--C-TTGGTTGA---TC---CCG
    Chimp	chr6	7000	7059	+	CGCTacccACTTTAgggTAGGATTAggaggTAT---AATTCTTaa-CTT---TGAGGGTCCAACCG

    """).lstrip()

    # Ensure it also works when the start offset isn't 0
    input.seek(0)
    block = Block.read_next_from_file(input)
    _, right_block = block.split(3)
    right_block.fill_in_first_sequence("CA")
    assert str(right_block) == dedent("""
    Mouse	chrX	53	55	+	CA------------
    Human	chr6	50	58	-	TGA---TC---CCG
    Chimp	chr6	7045	7059	+	TGAGGGTCCAACCG

    """).lstrip()

def test_merge_child_blocks_simple_overlap():
    """Run merge_child_blocks with two blocks that don't have perfect overlap."""
    input_a = StringIO(dedent("""
    foo	chr1	1	4	+	X-XX
    bar	chr1	51	55	-	ATGA

    """).lstrip())

    input_b = StringIO(dedent("""
    foo	chr1	1	6	+	XX-XXX
    baz	chr1	51	57	+	AGTTTA

    """).lstrip())

    chrom_sizes = {'chr1': 9, 'chr2': 3}

    output = StringIO()
    merge_child_blocks('foo', chrom_sizes, ['bar', 'baz'], [input_a, input_b], output)
    assert output.getvalue() == dedent("""
    foo	chr1	0	1	+	X

    foo	chr1	1	4	+	X-X-X
    bar	chr1	51	55	-	ATG-A
    baz	chr1	51	55	+	A-GTT

    foo	chr1	4	6	+	XX
    baz	chr1	55	57	+	TA

    foo	chr1	6	9	+	XXX

    foo	chr2	0	3	+	XXX

    """).lstrip()

def test_merge_child_blocks_ref_dups():
    """Make sure merge_child_blocks merges dups within the same blocks file."""
    input_a = StringIO(dedent("""
    foo	chr1	0	6	+	XX-XXXX
    bar	chr1	50	55	-	CATGA--

    foo	chr1	1	5	+	XXX-X
    bar	chr2	50	55	+	ATGCC

    """).lstrip())

    input_b = StringIO(dedent("""
    foo	chr1	0	6	+	XXX-XXX
    baz	chr1	50	57	+	CAGTTTA

    """).lstrip())

    chrom_sizes = {'chr1': 6}

    output = StringIO()
    merge_child_blocks('foo', chrom_sizes, ['bar', 'baz'], [input_a, input_b], output)

    # Columns:
    # 0 XC-C
    # 1 XAAA
    #   -T--
    # 2 XGTG
    #   ---T
    # 3 XAGT
    #   --C-
    # 4 X-CT
    # 5 X--A

    assert output.getvalue() == dedent("""
    foo	chr1	0	1	+	X
    bar	chr1	54	55	-	C
    baz	chr1	50	51	+	C

    foo	chr1	1	5	+	X-X-X-X
    bar	chr1	50	54	-	ATG-A--
    bar	chr2	50	55	+	A-T-GCC
    baz	chr1	51	56	+	A-GTT-T

    foo	chr1	5	6	+	X
    baz	chr1	56	57	+	A

    """).lstrip()

def test_merge_blocks_split_point_on_different_chrom():
    """Check that block ending in the same position on a different chrom does not get merged.

    This is a regression test for a bug which only triggered on a polytomy."""
    input_a = StringIO(dedent("""
    Anc01	Anc01refChr237	2	9	+	XXXXXXX
    ACERANA	KZ288318.1	138	145	+	AACtggt
    ACERANA	KZ288398.1	211122	211129	-	AACTGGt

    """).lstrip())

    input_b = StringIO(dedent("""
    Anc01	Anc01refChr237	2	9	+	XXXXXXX
    AMELLIFERA	CM009938.2	597321	597328	+	AACTGGT
    AMELLIFERA	CM009938.2	579437	579444	-	AACTGGT

    """).lstrip())

    input_c = StringIO(dedent("""
    Anc01	Anc01refChr238	0	9	+	XXXXXXXXX
    Anc00	Anc00refChr699	4813	4822	+	gatgattat

    """).lstrip())

    chrom_sizes = {'Anc01refChr238': 9, 'Anc01refChr237': 9}

    output = StringIO()
    merge_child_blocks('Anc01', chrom_sizes, ['foo', 'bar', 'baz'], [input_a, input_b, input_c], output)

    assert output.getvalue() == dedent("""
    Anc01	Anc01refChr237	0	2	+	XX

    Anc01	Anc01refChr237	2	9	+	XXXXXXX
    ACERANA	KZ288318.1	138	145	+	AACtggt
    ACERANA	KZ288398.1	211122	211129	-	AACTGGt
    AMELLIFERA	CM009938.2	597321	597328	+	AACTGGT
    AMELLIFERA	CM009938.2	579437	579444	-	AACTGGT

    Anc01	Anc01refChr238	0	9	+	XXXXXXXXX
    Anc00	Anc00refChr699	4813	4822	+	gatgattat

    """).lstrip()

def test_block_merge():
    input_1 = StringIO(dedent("""
    Human	chr6	50	77	+	CGG---GA---TCAACCAAG--AAGAATTGCTATA---
    Chimp	chr6	7000	7032	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcct

    """).lstrip())

    input_2 = StringIO(dedent("""
    Human	chr6	55	90	+	---TCAACCAAG--AAGAA-TTGCTATAcctccTAATCCTA---
    Gorilla	chr6	7000	7037	-	CCCTCA---AAGttAAGAAATT---ATAcctccTAATCCTAccc

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2 = Block.read_next_from_file(input_2)

    merged_block = block_1.merge(block_2)

    assert str(merged_block) == dedent("""
    Human	chr6	50	90	+	CGG---GA------TCAACCAAG----AAGAA-TTGCTATA---cctccTAATCCTA---
    Chimp	chr6	7000	7032	-	CGGTTGGACCC---TCA---AAGtt--AAGAA-TT---ATAcct----------------
    Gorilla	chr6	7000	7037	-	-----------CCCTCA---AAG--ttAAGAAATT---ATA---cctccTAATCCTAccc

    """).lstrip()

    # Try the other way around -- the results should be similar, except the
    # order of doubly-inserted regions should change.
    merged_block = block_2.merge(block_1)

    assert str(merged_block) == dedent("""
    Human	chr6	50	90	+	CGG---GA------TCAACCAAG----AAGAA-TTGCTATA---cctccTAATCCTA---
    Gorilla	chr6	7000	7037	-	--------CCC---TCA---AAGtt--AAGAAATT---ATA---cctccTAATCCTAccc
    Chimp	chr6	7000	7032	-	CGGTTGGA---CCCTCA---AAG--ttAAGAA-TT---ATAcct----------------

    """).lstrip()

def test_block_merge_adjacent():
    """Make sure Block.merge() works when the blocks don't overlap, but are adjacent."""
    input_1 = StringIO(dedent("""
    Human	chr6	50	77	+	CGG---GA---TCAACCAAG--AAGAATTGCTATA---
    Chimp	chr6	7000	7032	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcct

    """).lstrip())

    input_2 = StringIO(dedent("""
    Human	chr6	77	80	+	ATG

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2 = Block.read_next_from_file(input_2)

    merged_block = block_1.merge(block_2)

    assert str(merged_block) == dedent("""
    Human	chr6	50	80	+	CGG---GA---TCAACCAAG--AAGAATTGCTATA---ATG
    Chimp	chr6	7000	7032	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcct---

    """).lstrip()

def test_block_overlaps():
    input_1 = StringIO(dedent("""
    Human	chr6	50	77	+	CGG---GA---TCAACCAAG--AAGAATTGCTATA---
    Chimp	chr6	7000	7032	-	CGGTTGGACCCTCA---AAGttAAGAATT---ATAcct

    """).lstrip())

    input_2 = StringIO(dedent("""
    Human	chr6	76	79	+	XXX

    """).lstrip())

    input_3 = StringIO(dedent("""
    Human	chr6	78	7000	+	XXX
    Gorilla	chr6	7000	7032	-	XXX

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2 = Block.read_next_from_file(input_2)
    block_3 = Block.read_next_from_file(input_3)

    # Self should always overlap self, duh
    assert block_1.overlaps(block_1)

    # 2<->1 should overlap
    assert block_1.overlaps(block_2)
    assert block_2.overlaps(block_1)

    # Same with 2<->3
    assert block_2.overlaps(block_3)
    assert block_3.overlaps(block_2)

    # But 1<->3 shouldn't overlap
    assert not block_3.overlaps(block_1)
    assert not block_1.overlaps(block_3)

def test_block_merge_bug_with_offsets():
    """
    Test that Block.merge() still works even with BlockLines that have non-zero offset.
    """
    input_1 = StringIO(dedent("""
    mr	mrrefChr0	0	176	+	-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simMouse_chr6	simMouse.chr6	0	177	+	TTTTTCAGTTGCAATACCCAACCGGGAGAAACTTTCAGTGAGCACACCTCAGGTTCCTATATCAAGCAGGCAGTCTTGCATAGCAAATGGTCTCTGGTAGACGGTGCACTCAATCTATGTGAGGTATAGAAAATAAAGGACTACACACATCTCATCAAGTATCCCGTCATATTTGTG

    """).lstrip())

    input_2 = StringIO(dedent("""
    mr	mrrefChr0	146	176	+	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	144	174	+	CATCTCATCAAGTATCCTGTCGTATTTGTG

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2= Block.read_next_from_file(input_2)

    # Split the block to ensure the right block has non-zero start offsets
    _, split = block_1.split(144)

    # Check if this hits an assertion
    split.merge(block_2)

def test_merged_dup_stream_blocks():
    """Check merged_dup_stream behavior on Blocks."""
    input_1 = StringIO(dedent("""
    mr	mrrefChr0	12451	12577	+	----XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	12764	12894	+	ATTGATCTTATAGGAGCACTCCCTCTAGTCAGCAATATAATTAGATAGGACTGACACATTTGTAATCTTGTGCCTAATAGCAAAGGTGTGTGACGTACATCTAGATTGTGTAATGACCTGTAGTGTCAAC

    mr	mrrefChr0	12455	12456	+	X
    simRat_chr6	simRat.chr6	547138	547139	+	t

    mr	mrrefChr0	12458	12459	+	X
    simRat_chr6	simRat.chr6	547139	547140	+	a

    mr	mrrefChr0	12461	12464	+	XXX
    simRat_chr6	simRat.chr6	547140	547143	+	agt

    mr	mrrefChr0	12465	12503	+	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	547143	547181	+	ttgccTCTAGTCAGCAATCTAATTAGATAAGACTGACA

    mr	mrrefChr0	12503	12577	+	-------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	547181	547262	+	TACCGAACACTTGTAATCTTGTGCTTTTTAACAAAGGTGTATGTCATACATATAAATTGTGTGTTGACCTGTAGTGTCCAC

    """).lstrip())
    stream = merged_dup_stream(input_1, Block.read_next_from_file)

    block_1 = next(stream)
    assert str(block_1) == dedent("""
    mr	mrrefChr0	12451	12455	+	----XXXX
    simRat_chr6	simRat.chr6	12764	12772	+	ATTGATCT

    """).lstrip()

    block_2 = next(stream)
    assert str(block_2) == dedent("""
    mr	mrrefChr0	12455	12456	+	X
    simRat_chr6	simRat.chr6	12772	12773	+	T
    simRat_chr6	simRat.chr6	547138	547139	+	t

    """).lstrip()

    block_3 = next(stream)
    assert str(block_3) == dedent("""
    mr	mrrefChr0	12456	12458	+	XX
    simRat_chr6	simRat.chr6	12773	12775	+	AT

    """).lstrip()

    block_4 = next(stream)
    assert str(block_4) == dedent("""
    mr	mrrefChr0	12458	12459	+	X
    simRat_chr6	simRat.chr6	12775	12776	+	A
    simRat_chr6	simRat.chr6	547139	547140	+	a

    """).lstrip()

    block_5 = next(stream)
    assert str(block_5) == dedent("""
    mr	mrrefChr0	12459	12461	+	XX
    simRat_chr6	simRat.chr6	12776	12778	+	GG

    """).lstrip()

    block_6 = next(stream)
    assert str(block_6) == dedent("""
    mr	mrrefChr0	12461	12464	+	XXX
    simRat_chr6	simRat.chr6	12778	12781	+	AGC
    simRat_chr6	simRat.chr6	547140	547143	+	agt

    """).lstrip()

    block_7 = next(stream)
    assert str(block_7) == dedent("""
    mr	mrrefChr0	12464	12465	+	X
    simRat_chr6	simRat.chr6	12781	12782	+	A

    """).lstrip()

    block_8 = next(stream)
    assert str(block_8) == dedent("""
    mr	mrrefChr0	12465	12503	+	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	12782	12820	+	CTCCCTCTAGTCAGCAATATAATTAGATAGGACTGACA
    simRat_chr6	simRat.chr6	547143	547181	+	ttgccTCTAGTCAGCAATCTAATTAGATAAGACTGACA

    """).lstrip()

    block_10 = next(stream)
    assert str(block_10) == dedent("""
    mr	mrrefChr0	12503	12577	+	-------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simRat_chr6	simRat.chr6	12820	12894	+	-------CATTTGTAATCTTGTGCCTAATAGCAAAGGTGTGTGACGTACATCTAGATTGTGTAATGACCTGTAGTGTCAAC
    simRat_chr6	simRat.chr6	547181	547262	+	TACCGAACACTTGTAATCTTGTGCTTTTTAACAAAGGTGTATGTCATACATATAAATTGTGTGTTGACCTGTAGTGTCCAC

    """).lstrip()

    with pytest.raises(StopIteration):
        next(stream)

def test_merged_dup_stream_alignments():
    """Check merged_dup_stream behavior on AlignmentGroups."""
    alignments = StringIO(dedent("""
    chr5	20	25	Anc1c5	50	55	-
    chr6	0	2	Anc1c1	60	62	+
    chr6	0	3	Anc1c2	30	33	+
    chr6	4	11	Anc1c3	90	97	+
    chr6	5	10	Anc1c4	100	105	+
    """).lstrip())

    stream = merged_dup_stream(alignments, read_next_alignment)

    alignment_group_1 = next(stream)
    assert alignment_group_1.alignments == [Alignment('chr5', 20, 25, 'Anc1c5', 50, 55, '-')]

    alignment_group_2 = next(stream)
    assert alignment_group_2.alignments == [Alignment('chr6', 0, 2, 'Anc1c1', 60, 62, '+'),
                                            Alignment('chr6', 0, 2, 'Anc1c2', 30, 32, '+')]

    alignment_group_3 = next(stream)
    assert alignment_group_3.alignments == [Alignment('chr6', 2, 3, 'Anc1c2', 32, 33, '+')]

    alignment_group_4 = next(stream)
    assert alignment_group_4.alignments == [Alignment('chr6', 4, 5, 'Anc1c3', 90, 91, '+')]

    alignment_group_5 = next(stream)
    assert alignment_group_5.alignments == [Alignment('chr6', 5, 10, 'Anc1c3', 91, 96, '+'),
                                            Alignment('chr6', 5, 10, 'Anc1c4', 100, 105, '+')]

    alignment_group_6 = next(stream)
    assert alignment_group_6.alignments == [Alignment('chr6', 10, 11, 'Anc1c3', 96, 97, '+')]

    with pytest.raises(StopIteration):
        next(stream)

def test_block_compatible_with():
    input_1 = StringIO(dedent("""
    Human	chr6	0	2	+	A--T
    Chimp	chr6	8	12	+	ACCT

    """).lstrip())

    input_2 = StringIO(dedent("""
    Human	chr6	2	6	+	CCCC
    Chimp	chr6	12	16	+	GGGG

    """).lstrip())

    input_3 = StringIO(dedent("""
    Human	chr6	2	6	+	CCCC
    Chimp	chr6	12	16	-	GGGG

    """).lstrip())

    input_4 = StringIO(dedent("""
    Human	chr6	2	6	+	CCCC
    Gorilla	chr19	0	4	+	GGGG

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2 = Block.read_next_from_file(input_2)
    block_3 = Block.read_next_from_file(input_3)
    block_4 = Block.read_next_from_file(input_4)

    assert block_1.compatible_with(block_2)
    assert not block_2.compatible_with(block_1)

    assert not block_1.compatible_with(block_3)
    assert not block_2.compatible_with(block_3)

    assert block_1.compatible_with(block_4)

def test_block_concatenate():
    input_1 = StringIO(dedent("""
    Human	chr6	0	2	+	A------T
    Chimp	chr6	8	12	+	AC----CT
    Chimp	chr6	22	30	-	ACTGAGCT

    """).lstrip())

    input_2 = StringIO(dedent("""
    Human	chr6	2	6	+	C-CCC
    Chimp	chr6	18	22	-	GTGG-
    Chimp	chr6	12	16	+	GTG-G

    """).lstrip())

    block_1 = Block.read_next_from_file(input_1)
    block_2 = Block.read_next_from_file(input_2)

    block_1.concatenate(block_2)

    assert str(block_1) == dedent("""
    Human	chr6	0	6	+	A------TC-CCC
    Chimp	chr6	8	16	+	AC----CTGTG-G
    Chimp	chr6	18	30	-	ACTGAGCTGTGG-

    """).lstrip()

    # Test missing line partners on both sides (should be gapped out)
    input_3 = StringIO(dedent("""
    Human	chr6	0	2	+	A------T
    Chimp	chr6	8	12	+	AC----CT

    """).lstrip())

    input_4 = StringIO(dedent("""
    Human	chr6	2	6	+	C-CCC
    Gorilla	chr6	100	105	-	CCCCC

    """).lstrip())

    block_3 = Block.read_next_from_file(input_3)
    block_4 = Block.read_next_from_file(input_4)

    block_3.concatenate(block_4)

    assert str(block_3) == dedent("""
    Human	chr6	0	6	+	A------TC-CCC
    Chimp	chr6	8	12	+	AC----CT-----
    Gorilla	chr6	100	105	-	--------CCCCC

    """).lstrip()

def test_maximize_block_length():
    input = StringIO(dedent("""
    Human	chr6	0	2	+	X------X
    Chimp	chr6	8	12	+	AC----CT

    Human	chr6	2	6	+	X-XXX
    Gorilla	chr6	100	105	-	CCCCC

    Human	chr6	6	8	+	XX

    Human	chr6	8	12	+	XXXX
    Chimp	chr12	0	4	+	ACGT

    """).lstrip())

    ref_sequence = {'chr6': 'ACGTACGTACGT'}
    output = StringIO()
    maximize_block_length(input, ref_sequence, output)

    assert output.getvalue() == dedent("""
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    Human	chr6	8	12	+	ACGT
    Chimp	chr12	0	4	+	ACGT

    """).lstrip()

def test_block_to_maf():
    input = StringIO(dedent("""
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    """).lstrip())
    chrom_sizes = {'Human': {'chr6': 12},
                   'Chimp': {'chr6': 18},
                   'Gorilla': {'chr6': 200}}

    block = Block.read_next_from_file(input)
    maf_block = block.to_maf(chrom_sizes)

    assert maf_block == dedent("""
    a
    s Human.chr6 0 8 + 12 A------CG-TACGT
    s Chimp.chr6 8 4 + 18 AC----CT-------
    s Gorilla.chr6 95 5 - 200 --------CCCCC--

    """).lstrip()

def test_flip_alignment():
    input = StringIO(dedent("""
    simHuman.chr6	601632	601633	Anc1refChr1	332920	332921	+
    simHuman.chr6	601633	601639	Anc1refChr1	332921	332927	+
    simHuman.chr6	601639	601664	Anc1refChr1	332927	332952	-
    simHuman.chr6	601664	601670	Anc1refChr1	332952	332958	-
    simHuman.chr6	601670	601674	Anc1refChr1	332958	332962	+
    simHuman.chr6	601674	601684	Anc1refChr1	332963	332973	+
    simHuman.chr6	601685	601696	Anc1refChr1	332973	332984	+
    simHuman.chr6	601696	601703	Anc1refChr1	332984	332991	+
    simHuman.chr6	601703	601714	Anc1refChr1	333002	333013	+
    simHuman.chr6	601715	601827	Anc1refChr1	333013	333125	+
    """))
    output = StringIO()

    flip_alignment(input, output)

    assert output.getvalue() == dedent("""
    Anc1refChr1	332920	332921	simHuman.chr6	601632	601633	+
    Anc1refChr1	332921	332927	simHuman.chr6	601633	601639	+
    Anc1refChr1	332927	332952	simHuman.chr6	601639	601664	-
    Anc1refChr1	332952	332958	simHuman.chr6	601664	601670	-
    Anc1refChr1	332958	332962	simHuman.chr6	601670	601674	+
    Anc1refChr1	332963	332973	simHuman.chr6	601674	601684	+
    Anc1refChr1	332973	332984	simHuman.chr6	601685	601696	+
    Anc1refChr1	332984	332991	simHuman.chr6	601696	601703	+
    Anc1refChr1	333002	333013	simHuman.chr6	601703	601714	+
    Anc1refChr1	333013	333125	simHuman.chr6	601715	601827	+
    """).lstrip()

def test_lift_blocks():
    # simple
    alignments = StringIO(dedent("""
    chr6	0	2	Anc1refChr1	60	62	+
    chr6	5	7	Anc1refChr4	100	102	+
    """).lstrip())
    blocks = StringIO(dedent("""
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    Human	chr6	8	12	+	ACGT
    Chimp	chr12	0	4	+	ACGT

    """).lstrip())
    other_genome = 'Anc1'
    output = StringIO()

    lift_blocks(alignments, blocks, other_genome, output)

    assert output.getvalue() == dedent("""
    Anc1	Anc1refChr1	60	62	+	X------X----
    Human	chr6	0	5	+	A------CG-TA
    Chimp	chr6	8	12	+	AC----CT----
    Gorilla	chr6	101	105	-	--------CCCC

    Anc1	Anc1refChr4	100	102	+	XX-
    Human	chr6	5	8	+	CGT
    Gorilla	chr6	100	101	-	C--

    """).lstrip()

def test_alignment_group_split():
    alignments = [
        Alignment('chr1', 20, 60, 'otherChr1', 0, 40, '+'),
        Alignment('chr1', 20, 60, 'otherChr2', 20, 60, '-'),
    ]
    alignment_group = AlignmentGroup(alignments)
    left_group, right_group = alignment_group.split(20)
    assert left_group.chrom == 'chr1'
    assert left_group.start == 20
    assert left_group.end == 40
    assert left_group.alignments == [
        Alignment('chr1', 20, 40, 'otherChr1', 0, 20, '+'),
        Alignment('chr1', 20, 40, 'otherChr2', 40, 60, '-'),
    ]
    assert right_group.chrom == 'chr1'
    assert right_group.start == 40
    assert right_group.end == 60
    assert right_group.alignments == [
        Alignment('chr1', 40, 60, 'otherChr1', 20, 40, '+'),
        Alignment('chr1', 40, 60, 'otherChr2', 20, 40, '-'),
    ]


def test_lift_blocks_dups():
    alignments = StringIO(dedent("""
    chr6	0	2	Anc1c1	60	62	+
    chr6	0	3	Anc1c2	30	33	+
    chr6	4	11	Anc1c3	90	97	+
    chr6	5	10	Anc1c4	100	105	+
    """).lstrip())
    blocks = StringIO(dedent("""
    Human	chr6	0	8	+	A------CG-TACGT
    Chimp	chr6	8	12	+	AC----CT-------
    Gorilla	chr6	100	105	-	--------CCCCC--

    Human	chr6	8	12	+	ACGT
    Chimp	chr12	0	4	+	ACGT

    """).lstrip())
    other_genome = 'Anc1'
    output = StringIO()

    lift_blocks(alignments, blocks, other_genome, output)

    assert output.getvalue() == dedent("""
    Anc1	Anc1c1	60	62	+	X------X
    Human	chr6	0	2	+	A------C
    Chimp	chr6	8	12	+	AC----CT

    Anc1	Anc1c2	30	32	+	X------X
    Human	chr6	0	2	+	A------C
    Chimp	chr6	8	12	+	AC----CT

    Anc1	Anc1c2	32	33	+	X--
    Human	chr6	2	4	+	G-T
    Gorilla	chr6	102	105	-	CCC

    Anc1	Anc1c3	90	91	+	X
    Human	chr6	4	5	+	A
    Gorilla	chr6	101	102	-	C

    Anc1	Anc1c3	91	94	+	XXX
    Human	chr6	5	8	+	CGT
    Gorilla	chr6	100	101	-	C--

    Anc1	Anc1c4	100	103	+	XXX
    Human	chr6	5	8	+	CGT
    Gorilla	chr6	100	101	-	C--

    Anc1	Anc1c3	94	96	+	XX
    Human	chr6	8	10	+	AC
    Chimp	chr12	0	2	+	AC

    Anc1	Anc1c4	103	105	+	XX
    Human	chr6	8	10	+	AC
    Chimp	chr12	0	2	+	AC

    Anc1	Anc1c3	96	97	+	X-
    Human	chr6	10	12	+	GT
    Chimp	chr12	2	4	+	GT

    """).lstrip()

def test_merged_dup_stream_same_starts():
    """
    Test merged_dup_stream against a rare bug that crops up with a pattern of dups.
    """
    block_file = StringIO(dedent("""
Anc2	Anc2refChr0	11463	11531	+	---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
simCow_chr6	simCow.chr6	129132	130571	+	tgaagttaaagggtggtggaaggactaccaaagcaatgggaaaggcaaaaagaaaacgggtcggatcgcaatcccaattaccaagaagttggactatgatgttaagtcagcaagagtcaaaaggaaaagagtcaaggaagggcattcagtgatggtcagaggctctattcaacctgagagattaaccgtagagaatgtacacaaacccaatacaggtactaagactgagcatgcaagggatctgtgtcttgactcatataccgtagtagtggagcatgcgaacacaccggcatcgaccctggataagtcaacacaaagaaaagttaacacttacacacacgagctgaactcagttttgcatcaagccaacttgatcgacatctaccaaactttacatcctattccggctggttactcatttttttctgatccataccattcatacttaaagatagattacatcgtcgtgaaagcactcgtgcataaatgcaaagactccgatatagtagctaattgcctcgctgatcagagctcaatgaagttgagtctggagataagtgctctcacccagaatcagttgactccaaattgctcaaattggcatccgaataatgtatcattccttctaactgagtattgtgtccacaatgagatgaagacagaaattaagacgtttcttgacaaaaactaaaattgggatctcacctatggaaatacttgggagattggggcgttaatcatgctaaacaacacaggtaaaactgggcaaggaaggttgaggatggtgatccttaagtcaaaaataaaagaggttcattagcaagagcttacgcacgccagggcggctaggcagcaacaaatgactaaaataaaatttgagctgaagagtgtaaagaccgaagacacattgcccaaggtcaatgaatccaagagctggtttcttgacaaaaacaaaatcgacaggtctttgctgagcgtgcagaatgtagaaaaaatcatatcgacgtgttgaggaatcatcaggggtacataacaaccgacccagtggacattcaaataggcattaagcagtatcaccaacctctttgcaaaagtttaacttagttgacattgatattttcgattcgtttacccttctgcgactcaacaaagaaggaatagagtctctaaatagatcagtaatcggccccgaaataaagcatttattaacagtccatcgaagaagttacctggtttgccactgatattttacaagcagtacaaagaaaagttggtctcatttttttcctcaagctacaatcaattgagaatgaaggaatactaccgaactattttcacgagtcgtgtagtgtatcgatatcgaaatcaagcgaagccatgatgaaaagacaaaaggtgacatggcccactttttagacatcacacgccggaaaactaaa

Anc2	Anc2refChr0	11463	12315	+	----XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
simCow_chr6	simCow.chr6	11254	12110	+	CAGATaagccatgatgaaaagacagaaggcgacgcggcccgttttttcgacatcacacgatggaaaactaaacagaatagtagccaaccgaaactatcaacgaatcggaaaacttatacattatgatcacttaggtttcaccccagggatgccaggatggcttaggatacgcaaacatataaatgtaatttagcctgttaaccaaactaatcagaaaaattacatgatactatcattagacagggatgatggaccgataaaacaaactgacccagagtacatgcttacaatcctgaatgaattggggattgaaggcaattacctcaatgttgtgagggccatttatgataaaccatcatccaatctcatagtcaacgggcgtaagcttgaagtctttcccttggaaggcaccaggaggcggcaaggttgctctatctctccactcttgttcaacatagtgctagaattgttggcccgcactagccgtcatgagagagacgagcggggtatacaattggatgaagaggaaacaaagctagccatcttcgcggatgacatgatagtctacttggggaacccgatggcggcagcaggcaacctacttaaactaaaaagcaacttgagtgaactctctggttataagatcaacgtgtcaagtcccaaacattcttgtacaaaaacaaccagcagaccagtcatagtataagcgagtcacctccgaccaccgattctaaatgaatcaaattTTtaggaactcaactgacccaagacatcaaggacttgttgatccggatacggaacgtgctattaagataaataagggaaattactaataatccccacaatattaagtgtt

    Anc2	Anc2refChr0	11515	11531	+	XXXXXXXXXXXXXXXX
    simCow_chr6	simCow.chr6	464111	464127	-	atgatggaagtataaa

    """).lstrip())

    stream = merged_dup_stream(block_file, Block.read_next_from_file)

    block_1 = next(stream)
    assert str(block_1) == dedent("""
    Anc2	Anc2refChr0	11463	11515	+	-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    simCow_chr6	simCow.chr6	129132	130555	+	tgaagttaaagggtggtggaaggactaccaaagcaatgggaaaggcaaaaagaaaacgggtcggatcgcaatcccaattaccaagaagttggactatgatgttaagtcagcaagagtcaaaaggaaaagagtcaaggaagggcattcagtgatggtcagaggctctattcaacctgagagattaaccgtagagaatgtacacaaacccaatacaggtactaagactgagcatgcaagggatctgtgtcttgactcatataccgtagtagtggagcatgcgaacacaccggcatcgaccctggataagtcaacacaaagaaaagttaacacttacacacacgagctgaactcagttttgcatcaagccaacttgatcgacatctaccaaactttacatcctattccggctggttactcatttttttctgatccataccattcatacttaaagatagattacatcgtcgtgaaagcactcgtgcataaatgcaaagactccgatatagtagctaattgcctcgctgatcagagctcaatgaagttgagtctggagataagtgctctcacccagaatcagttgactccaaattgctcaaattggcatccgaataatgtatcattccttctaactgagtattgtgtccacaatgagatgaagacagaaattaagacgtttcttgacaaaaactaaaattgggatctcacctatggaaatacttgggagattggggcgttaatcatgctaaacaacacaggtaaaactgggcaaggaaggttgaggatggtgatccttaagtcaaaaataaaagaggttcattagcaagagcttacgcacgccagggcggctaggcagcaacaaatgactaaaataaaatttgagctgaagagtgtaaagaccgaagacacattgcccaaggtcaatgaatccaagagctggtttcttgacaaaaacaaaatcgacaggtctttgctgagcgtgcagaatgtagaaaaaatcatatcgacgtgttgaggaatcatcaggggtacataacaaccgacccagtggacattcaaataggcattaagcagtatcaccaacctctttgcaaaagtttaacttagttgacattgatattttcgattcgtttacccttctgcgactcaacaaagaaggaatagagtctctaaatagatcagtaatcggccccgaaataaagcatttattaacagtccatcgaagaagttacctggtttgccactgatattttacaagcagtacaaagaaaagttggtctcatttttttcctcaagctacaatcaattgagaatgaaggaatactaccgaactattttcacgagtcgtgtagtgtatcgatatcgaaatcaagc----gaagccatgatgaaaagacaaaaggtgacatggcccactttttagacatcac
    simCow_chr6	simCow.chr6	11254	11310	+	---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAGATaagccatgatgaaaagacagaaggcgacgcggcccgttttttcgacatcac

    """).lstrip()

    block_2 = next(stream)
    assert str(block_2) == dedent("""
    Anc2	Anc2refChr0	11515	11531	+	XXXXXXXXXXXXXXXX
    simCow_chr6	simCow.chr6	130555	130571	+	acgccggaaaactaaa
    simCow_chr6	simCow.chr6	11310	11326	+	acgatggaaaactaaa
    simCow_chr6	simCow.chr6	464111	464127	-	atgatggaagtataaa

    """).lstrip()

def test_lift_region():
    alignment_1 = Alignment(my_chrom='Anc1refChr1', my_start=91810, my_end=91881,
                          other_chrom='mrrefChr5', other_start=57149, other_end=57220,
                          other_strand='+')
    assert lift_region(alignment_1, 91810, 91881) == (57149, 57220)
    assert lift_region(alignment_1, 91845, 91870) == (57184, 57209)

    alignment_2 = Alignment(my_chrom='Anc1refChr0', my_start=8054, my_end=8060,
                          other_chrom='mrrefChr0', other_start=44139, other_end=44145,
                          other_strand='-')
    assert lift_region(alignment_2, 8054, 8060) == (44139, 44145)
    assert lift_region(alignment_2, 8056, 8059) == (44140, 44143)

def test_group_by_size():
    chrom_sizes = {'chr1': 5000, 'chr2': 10000, 'chrY': 1000}
    files_by_chroms = [{'chr1': NamedTemporaryFile().name, 'chr2': NamedTemporaryFile().name, 'chrY': NamedTemporaryFile().name}]
    for chrom, path in list(files_by_chroms[0].items()):
        with open(path, 'w') as f:
            f.write(chrom + "\n")
    chrom_sizes_subsets, files = group_by_size(files_by_chroms, chrom_sizes, lambda: NamedTemporaryFile().name, target_size=7000)
    assert chrom_sizes_subsets == [{'chr1': 5000, 'chrY': 1000}, {'chr2': 10000}]
    assert all([len(f) == 1 for f in files]), "Should have only one output file per grouping if given one input file"
    file_contents = [open(f[0]).read() for f in files]
    assert file_contents == ['chr1\nchrY\n', 'chr2\n']

def test_split_blocks_by_chrom():
    block_file = StringIO(dedent("""
    A	chr1	0	5	+	XXXXX
    B	chr1	0	5	+	ACGTG

    A	chr1	15	20	+	XXXXX
    B	chr1	20	25	+	ACGTG

    A	chr2	15	20	+	XXXXX
    B	chr1	20	25	+	ACGTG

    """).lstrip())

    paths_by_chrom = split_blocks_by_chrom(block_file, get_temp_path=lambda: NamedTemporaryFile().name)
    assert len(paths_by_chrom) == 2
    with open(paths_by_chrom['chr1']) as f:
        assert f.read() == dedent("""
        A	chr1	0	5	+	XXXXX
        B	chr1	0	5	+	ACGTG

        A	chr1	15	20	+	XXXXX
        B	chr1	20	25	+	ACGTG

        """).lstrip()
    with open(paths_by_chrom['chr2']) as f:
        assert f.read() == dedent("""
        A	chr2	15	20	+	XXXXX
        B	chr1	20	25	+	ACGTG

        """).lstrip()

def test_split_alignments_by_chrom():
    alignments_file = StringIO(dedent("""
    simHuman.chr6	601632	601633	Anc1refChr1	332920	332921	+
    simHuman.chr6	601633	601639	Anc1refChr1	332921	332927	+
    simHuman.chr6	601639	601664	Anc1refChr1	332927	332952	-
    simHuman.chrX	601664	601670	Anc1refChr1	332952	332958	-
    simHuman.chrX	601670	601674	Anc1refChr1	332958	332962	+
    simHuman.chrY	601674	601684	Anc1refChr1	332963	332973	+
    """).lstrip())
    paths_by_chrom = split_alignments_by_chrom(alignments_file, get_temp_path=lambda: NamedTemporaryFile().name)
    assert len(paths_by_chrom) == 3
    with open(paths_by_chrom['simHuman.chr6']) as f:
        assert f.read() == dedent("""
        simHuman.chr6	601632	601633	Anc1refChr1	332920	332921	+
        simHuman.chr6	601633	601639	Anc1refChr1	332921	332927	+
        simHuman.chr6	601639	601664	Anc1refChr1	332927	332952	-
        """).lstrip()
    with open(paths_by_chrom['simHuman.chrX']) as f:
        assert f.read() == dedent("""
        simHuman.chrX	601664	601670	Anc1refChr1	332952	332958	-
        simHuman.chrX	601670	601674	Anc1refChr1	332958	332962	+
        """).lstrip()
    with open(paths_by_chrom['simHuman.chrY']) as f:
        assert f.read() == dedent("""
        simHuman.chrY	601674	601684	Anc1refChr1	332963	332973	+
        """).lstrip()

if __name__ == '__main__':
    main()
