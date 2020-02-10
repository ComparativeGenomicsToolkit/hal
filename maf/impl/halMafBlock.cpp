/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafBlock.h"
#include <cstdlib>
#include <iostream>
#include <limits>

using namespace std;
using namespace hal;

const hal_index_t MafBlock::defaultMaxLength = 1000;

MafBlock::MafBlock(hal_index_t maxLength) : _maxLength(maxLength), _fullNames(false), _tree(NULL) {
    if (_maxLength <= 0) {
        _maxLength = numeric_limits<hal_index_t>::max();
    }
}

MafBlock::~MafBlock() {
    for (Entries::iterator i = _entries.begin(); i != _entries.end(); ++i) {
        delete i->second;
    }
    for (size_t j = 0; j < _stringBuffers.size(); ++j) {
        delete _stringBuffers[j];
    }
    if (_printTree && _tree != NULL) {
        stTree_destruct(_tree);
    }
}

void MafBlock::resetEntries() {
    _reference = NULL;
    _refIndex = NULL_INDEX;
    Entries::iterator i = _entries.begin();
    while (i != _entries.end()) {
        Entries::iterator next = i;
        ++next;
        MafBlockEntry *e = i->second;
        bool deleted = false;

        // every time we reset an entry, we check if was empty.
        // if it was, then we increase lastUsed, otherwise we reset it to
        // 0. this way, if an entry isn't used within 10 consecutive resets
        // we ditch it.  this is a tuning parameter to balance between not
        // creating and destroying entries every time there is a delete and
        // not keeping track of thousands of entries (fly scaffolds) which
        // don't get used but bog down all set operations in the flys.
        if (e->_start == NULL_INDEX) {
            if (e->_lastUsed > 10) {
                delete i->second;
                _entries.erase(i);
                deleted = true;
            } else {
                ++e->_lastUsed;
            }
        } else {
            e->_lastUsed = 0;
        }
        if (deleted == false) {
            assert(e->_start == NULL_INDEX || e->_length > 0);
            assert(e->_name == getName(i->first));
            // Rest block information but leave sequence information so we
            // can reuse it.
            e->_start = NULL_INDEX;
            e->_strand = '+';
            e->_length = 0;
            e->_sequence->clear();
        }
        i = next;
    }
}

void MafBlock::initEntry(MafBlockEntry *entry, const Sequence *sequence, DnaIteratorPtr dna, bool clearSequence) {
    string sequenceName = getName(sequence);
    if (entry->_name != sequenceName || sequence->getGenome() != entry->_genome) {
        // replace genearl sequence information
        entry->_name = sequenceName;
        entry->_genome = sequence->getGenome();
        entry->_srcLength = (hal_index_t)sequence->getSequenceLength();
    }
    if (dna != nullptr) {
        // update start position from the iterator
        entry->_start = dna->getArrayIndex() - sequence->getStartPosition();
        entry->_length = 0;
        entry->_strand = dna->getReversed() ? '-' : '+';
        if (dna->getReversed()) {
            entry->_start = entry->_srcLength - 1 - entry->_start;
        }
        if ((entry->_name  == "BTERRESTRIS.NC_015770.1") and (entry->_strand == '-') ) {  // FIXME hack
            if (false) {
                throw hal_exception("FIXME: Reference entry in MAF block created on negative strand: " + entry->_name);
            } else {
                cerr << "FIXME: Reference entry in MAF block created on negative strand: " << entry->_name << endl;
            }
        }
    } else {
        // no start position, so we wait till next time
        entry->_start = NULL_INDEX;
        entry->_length = 0;
        entry->_strand = '+';
    }
    if (clearSequence == true) {
        entry->_sequence->clear();
    }
    entry->_tree = NULL;
}

inline void MafBlock::updateEntry(MafBlockEntry *entry, const Sequence *sequence, DnaIteratorPtr dna) {
    if (dna != nullptr) {
        if (entry->_start == NULL_INDEX) {
            initEntry(entry, sequence, dna, false);
        }
        assert(entry->_name == getName(sequence));
        assert(entry->_strand == (dna->getReversed() ? '-' : '+'));
        assert(entry->_srcLength == (hal_index_t)sequence->getSequenceLength());

        ++entry->_length;

        assert(dna->getReversed() == true ||
               (hal_index_t)(dna->getArrayIndex() - sequence->getStartPosition()) ==
                   (hal_index_t)(entry->_start + entry->_length - 1));

        assert(dna->getReversed() == false ||
               (hal_index_t)(entry->_srcLength - 1 - (dna->getArrayIndex() - sequence->getStartPosition())) ==
                   (hal_index_t)(entry->_start + entry->_length - 1));

        entry->_sequence->append(dna->getBase());
    } else {
        entry->_sequence->append('-');
    }
}

// Puts the given node and its parents at the start of all their
// children lists. This has the effect of making the node first in a
// post-order traversal.
static void prioritizeNodeInTree(stTree *node) {
    stTree *parent = stTree_getParent(node);
    if (parent == NULL) {
        // Nothing to do.
        return;
    }

    // Swap the first node and this node.
    int64_t nodeIndex = -1;
    for (int64_t i = 0; i < stTree_getChildNumber(parent); i++) {
        if (stTree_getChild(parent, i) == node) {
            nodeIndex = i;
            break;
        }
    }
    assert(nodeIndex != -1);

    stTree *tmp = stTree_getChild(parent, 0);
    stTree_setChild(parent, 0, node);
    stTree_setChild(parent, nodeIndex, tmp);
    prioritizeNodeInTree(parent);
}

stTree *MafBlock::getTreeNode(SegmentIteratorPtr segIt, bool modifyEntries) {
    // Make sure the segment is sliced to only 1 base.
    assert(segIt->getStartPosition() == segIt->getEndPosition());
    stTree *ret = stTree_construct();
    const Genome *genome = segIt->getGenome();
    const Sequence *seq = genome->getSequenceBySite(segIt->getStartPosition());
    Entries::const_iterator entryIt = _entries.lower_bound(seq);
    if (entryIt != _entries.end() && entryIt->first == seq) {
        MafBlockEntry *entry = NULL;
        for (; entryIt != _entries.end() && entryIt->first == seq; entryIt++) {
            MafBlockEntry *curEntry = entryIt->second;
            hal_index_t curEntryPos = curEntry->_start + curEntry->_length;
            if (curEntry->_strand == '-') {
                curEntryPos = curEntry->_srcLength - 1 - curEntryPos;
            }
            if (curEntryPos == segIt->getStartPosition() - seq->getStartPosition() || curEntry->_start == NULL_INDEX) {
                entry = curEntry;
                break;
            }
        }
        assert(entry != NULL);
        stTree_setClientData(ret, entry);
        stTree_setLabel(ret, stString_copy(entry->_name.c_str()));
        if (modifyEntries) {
            entry->_tree = ret;
        }
    } else {
        // No entry for this sequence. Can happen if this is an ancestor
        // and we aren't including ancestral sequence.
        entryIt = _entries.begin();
        while (entryIt != _entries.end()) {
            entryIt++;
        }
        assert(genome->getNumChildren() != 0);
        stTree_setLabel(ret, stString_copy(segIt->getGenome()->getName().c_str()));
        stTree_setClientData(ret, NULL);
    }

    return ret;
}

// tree parameter represents node corresponding to the genome with
// bottom segment botIt
void MafBlock::buildTreeR(BottomSegmentIteratorPtr botIt, stTree *tree, bool modifyEntries) {
    const Genome *genome = botIt->getGenome();

    // attach a node and recurse for each of this segment's children
    // (and paralogous segments)
    for (hal_size_t i = 0; i < botIt->bseg()->getNumChildren(); i++) {
        if (botIt->bseg()->hasChild(i)) {
            const Genome *child = genome->getChild(i);
            TopSegmentIteratorPtr topIt = child->getTopSegmentIterator();
            topIt->toChild(botIt, i);
            stTree *canonicalParalog = getTreeNode(topIt, modifyEntries);
            stTree_setParent(canonicalParalog, tree);
            if (topIt->tseg()->hasParseDown()) {
                BottomSegmentIteratorPtr childBotIt = child->getBottomSegmentIterator();
                childBotIt->toParseDown(topIt);
                buildTreeR(childBotIt, canonicalParalog, modifyEntries);
            }
            // Traverse the paralogous segments cycle and add those segments as well
            assert(topIt->tseg()->isCanonicalParalog());
            if (topIt->tseg()->hasNextParalogy()) {
                topIt->toNextParalogy();
                while (!topIt->tseg()->isCanonicalParalog()) {
                    stTree *paralog = getTreeNode(topIt, modifyEntries);
                    stTree_setParent(paralog, tree);
                    if (topIt->tseg()->hasParseDown()) {
                        BottomSegmentIteratorPtr childBotIt = child->getBottomSegmentIterator();
                        childBotIt->toParseDown(topIt);
                        buildTreeR(childBotIt, paralog, modifyEntries);
                    }
                    topIt->toNextParalogy();
                }
            }
        }
    }
}

stTree *MafBlock::buildTree(ColumnIteratorPtr colIt, bool modifyEntries) {
    // Get any base from the column to begin building the tree
    const ColumnMap *colMap = colIt->getColumnMap();
    ColumnMap::const_iterator colMapIt = colMap->begin();
    const Sequence *sequence = NULL;
    hal_index_t index = NULL_INDEX;
    while (colMapIt != colMap->end()) {
        if (!colMapIt->second->empty()) {
            // Found a non-empty column map entry, just take the index and
            // sequence of the first base found
            sequence = colMapIt->first;
            index = colMapIt->second->at(0)->getArrayIndex();
            break;
        }
        colMapIt++;
    }
    assert(sequence != NULL && index != NULL_INDEX);
    const Genome *genome = sequence->getGenome();

    // Get the bottom segment that is the common ancestor of all entries
    TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIt;
    if (genome->getNumTopSegments() == 0) {
        // The reference is the root genome.
        botIt = genome->getBottomSegmentIterator();
        botIt->toSite(index);
    } else {
        // Keep heading up the tree until we hit the root segment.
        topIt->toSite(index);
        while (topIt->tseg()->hasParent()) {
            const Genome *parent = topIt->getGenome()->getParent();
            botIt = parent->getBottomSegmentIterator();
            botIt->toParent(topIt);
            if (parent->getParent() == NULL || !botIt->bseg()->hasParseUp()) {
                // Reached root genome
                break;
            }
            topIt = parent->getTopSegmentIterator();
            topIt->toParseUp(botIt);
        }
    }

    stTree *tree = NULL;
    if ((not topIt->tseg()->hasParent()) && (topIt->getGenome() == genome) && (genome->getNumBottomSegments() == 0)) {
        // Handle insertions in leaves. botIt doesn't point anywhere since
        // there are no bottom segments.
        tree = getTreeNode(topIt, modifyEntries);
    } else {
        tree = getTreeNode(botIt, modifyEntries);
        buildTreeR(botIt, tree, modifyEntries);
    }
    assert(tree != NULL);
    return tree;
}

// No DNA Iterator for this sequence.  We just give it an empty entry
void MafBlock::initBlockEmpty(ColumnMap::const_iterator& colMapIt) {
    const Sequence *sequence = colMapIt->first;
    Entries::iterator entryIt = _entries.lower_bound(sequence);
    if (entryIt == _entries.end() || entryIt->first != sequence) {
        MafBlockEntry *entry = new MafBlockEntry(_stringBuffers);
        initEntry(entry, sequence, DnaIteratorPtr());
        entryIt = _entries.insert(Entries::value_type(sequence, entry));
    } else {
        assert(entryIt->first == sequence);
        assert(entryIt->second->_name == getName(sequence));
        initEntry(entryIt->second, sequence, DnaIteratorPtr());
    }
}

// Have DNA Iterator for this sequence, so fill in the block
void MafBlock::initBlockFill(ColumnMap::const_iterator& colMapIt) {
    const Sequence *sequence = colMapIt->first;
    Entries::iterator entryIt = _entries.lower_bound(sequence);
    for (DNASet::const_iterator d = colMapIt->second->begin(); d != colMapIt->second->end(); ++d) {
        // search for colMapIt's sequence in _entries.
        // we conly call find() once.  afterwards we just move forward
        // in the map since they are both sorted by the same key.
        if (entryIt == _entries.begin()) {
            entryIt = _entries.lower_bound(sequence);
            if (entryIt == _entries.end() || entryIt->first != sequence) {
                entryIt = _entries.end();
            }
        } else {
            while ((entryIt->first != colMapIt->first) && (entryIt != _entries.end())) {
                ++entryIt;
            }
        }
        if (entryIt == _entries.end()) {
            MafBlockEntry *entry = new MafBlockEntry(_stringBuffers);
            initEntry(entry, sequence, *d);
            assert(entry->_name == getName(sequence));
            entryIt = _entries.insert(Entries::value_type(sequence, entry));
        } else {
            initEntry(entryIt->second, sequence, *d);
        }
        ++entryIt;
    }
}

// initialize the reference
void MafBlock::initBlockInitReference(ColumnIteratorPtr col) {
    const Sequence *referenceSequence = col->getReferenceSequence();
    Entries::iterator entryIt = _entries.lower_bound(referenceSequence);
    if (entryIt == _entries.end() || entryIt->first != referenceSequence) {
        entryIt = _entries.begin();
    }
    _reference = entryIt->second;
    if (entryIt->first == referenceSequence) {
        _refIndex = col->getReferenceSequencePosition();
    }
}

void MafBlock::initBlock(ColumnIteratorPtr col, bool fullNames, bool printTree) {
    if (printTree && _tree != NULL) {
        stTree_destruct(_tree);
    }
    resetEntries();
    _fullNames = fullNames;
    _printTree = printTree;
    const ColumnMap *colMap = col->getColumnMap();
    static int cnt_FIXME = 0;

    for (ColumnMap::const_iterator colMapIt = colMap->begin(); colMapIt != colMap->end(); ++colMapIt) {
        const Sequence *sequence = colMapIt->first; // FIXME
        if (sequence->getName() == "NC_015770.1") { // FIXME
            cnt_FIXME++;
            cerr << "initBlock on " << sequence->getName() << " " << cnt_FIXME << endl;
        }

        if (colMapIt->second->empty()) {
            initBlockEmpty(colMapIt);
        } else {
            initBlockFill(colMapIt);
        }
    }

    if (_reference == NULL) {
        initBlockInitReference(col);
    }

    if (_printTree) {
        _tree = buildTree(col, true);
    }
}

void MafBlock::appendColumn(ColumnIteratorPtr col) {
    const ColumnMap *colMap = col->getColumnMap();
    Entries::iterator e = _entries.begin();
    ColumnMap::const_iterator c = colMap->begin();
    DNASet::const_iterator d;
    const Sequence *sequence;

    for (; c != colMap->end(); ++c) {
        sequence = c->first;
        for (d = c->second->begin(); d != c->second->end(); ++d) {
            while (e->first != sequence && e != _entries.end()) {
                updateEntry(e->second, NULL, DnaIteratorPtr());
                ++e;
            }
            assert(e != _entries.end());
            assert(e->first == sequence);
            assert(e->second->_name == getName(sequence));
            updateEntry(e->second, sequence, *d);
            ++e;
        }
    }

    for (; e != _entries.end(); ++e) {
        updateEntry(e->second, NULL, DnaIteratorPtr());
    }
}

// Q: When can we append a column?
// A: When for every sequence already in the column, the new column
//    has either a gap or a contigugous base.  The new column also has
//    no new sequences.
bool MafBlock::canAppendColumn(ColumnIteratorPtr col) {
    const ColumnMap *colMap = col->getColumnMap();
    Entries::iterator e = _entries.begin();

    for (ColumnMap::const_iterator c = colMap->begin(); c != colMap->end(); ++c) {
        const Sequence *sequence = c->first;
        hal_index_t sequenceStart = sequence->getStartPosition();

        for (DNASet::const_iterator d = c->second->begin(); d != c->second->end(); ++d) {
            while (e->first != sequence && e != _entries.end()) {
                ++e;
            }
            if (e == _entries.end()) {
                return false;
            } else {
                MafBlockEntry *entry = e->second;
                assert(e->first == sequence);
                assert(entry->_name == getName(sequence) && entry->_genome == sequence->getGenome());
                if (entry->_start != NULL_INDEX) {
                    if (entry->_length >= _maxLength ||
                        (entry->_length > 0 && (entry->_strand == '-') != (*d)->getReversed())) {
                        return false;
                    }
                    hal_index_t pos = (*d)->getArrayIndex() - sequenceStart;
                    if ((*d)->getReversed() == true) {
                        // position on reverse strand relative to end of sequence
                        pos = entry->_srcLength - 1 - pos;
                    }
                    if (pos - entry->_start != entry->_length) {
                        return false;
                    }
                }
                ++e;
            }
        }
    }
    if (_printTree) {
        stTree *tree = buildTree(col, false);
        bool ret = stTree_equals(tree, _tree);
        stTree_destruct(tree);
        return ret;
    }
    return true;
}

ostream &hal::operator<<(ostream &os, const MafBlockEntry &mafBlockEntry) {
    os << "s\t" << mafBlockEntry._name << '\t' << mafBlockEntry._start << '\t' << mafBlockEntry._length << '\t'
       << mafBlockEntry._strand << '\t' << mafBlockEntry._srcLength << '\t' << mafBlockEntry._sequence->str() << '\n';
    return os;
}

istream &hal::operator>>(istream &is, MafBlockEntry &mafBlockEntry) {
    string buffer;
    is >> mafBlockEntry._name >> mafBlockEntry._start >> mafBlockEntry._length >> mafBlockEntry._strand >>
        mafBlockEntry._srcLength >> buffer;
    assert(mafBlockEntry._strand == '+' || mafBlockEntry._strand == '-');
    // don't think this fucntion is used so don't worry about
    // this crap too much.
    mafBlockEntry._sequence->clear();
    for (size_t i = 0; i < buffer.length(); ++i) {
        mafBlockEntry._sequence->append(buffer[i]);
    }
    return is;
}

static void printTreeEntries(stTree *tree, ostream &os) {
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        printTreeEntries(child, os);
    }
    MafBlockEntry *entry = (MafBlockEntry *)stTree_getClientData(tree);
    if (entry != NULL) {
        // The entry can be null if --noAncestors is enabled.
        os << *entry;
    }
}

void MafBlock::checkRefStrand() const {
    if ((_reference != NULL) and (_reference->_strand == '-')) {
        throw hal_exception("Reference entry in MAF block created on negative strand: " + _reference->_name);
    }
}

ostream &MafBlock::printBlockWithTree(ostream &os) const {
    // Sort tree so that the reference comes first.
    prioritizeNodeInTree(_reference->_tree);

    // Print tree as a block comment.
    char *treeString = stTree_getNewickTreeString(_tree);
    os << "a tree=\"" << treeString << "\"\n";

    // Print entries in post order.
    printTreeEntries(_tree, os);
    checkRefStrand(); // after block is written
    return os;
}

// todo: fast way of reference first.
ostream &MafBlock::printBlock(ostream &os) const {
    os << "a\n";

    assert(_reference != NULL);
    if (_reference->_start == NULL_INDEX) {
        if (_refIndex != NULL_INDEX) {
            _reference->_start = _refIndex;
            os << *_reference;
            _reference->_start = NULL_INDEX;
        }
    } else {
        os << *_reference;
    }

    for (Entries::const_iterator e = _entries.begin(); e != _entries.end(); ++e) {
        if ((e->second->_start != NULL_INDEX) && (e->second != _reference)) {
            os << *e->second;
        }
    }
    checkRefStrand(); // after block is written so it can be looked at.
    return os;
}

ostream &hal::operator<<(ostream &os, const MafBlock &mafBlock) {
    if (mafBlock._printTree) {
        return mafBlock.printBlockWithTree(os);
    } else {
        return mafBlock.printBlock(os);
    }
}

istream &hal::operator>>(istream &is, MafBlock &mafBlock) {
    return is;
}
