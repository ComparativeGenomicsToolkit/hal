#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

'''
Objects & functions to parse, manipulate, draw, or get info from, trees
'''

from Bio import Phylo
from Bio.Phylo.Newick import Clade
import os, copy, sys
from sonLib.bioio import system  

def isBinaryTree(tree):
    for clade in tree.get_nonterminals():
        children = clade.clades
        if len(children) != 2:
            #sys.stderr.write("NODE %s has %d children!!\n" %(clade.name, len(children)))
            #print children
            #sys.stderr.write("Children: %s\n" %(",".join([c.name for c in children])))
            #sys.exit(1)
            return False
    return True

def getLeftRight(children):
    #Left child has less leaves than right child
    count0 = children[0].count_terminals()
    count1 = children[1].count_terminals()
    if count0 <= count1:
        return children[0], children[1]
    else:
        return children[1], children[0] 

def inorder(tree, reverse=False):
    #Return a middle-sorted (left, root, right) list of node names
    if tree.is_terminal():
        if tree.name:
            return [tree.name]
        else:
            return []
    rootName = []
    if tree.name:
        rootName = [tree.name]
    children = tree.clades
    if len(children) < 2:
        return inorder(children[0]) + rootName
    left, right = getLeftRight(children)
    leftNames = inorder(left)
    rightNames = inorder(right)
    if not reverse:
        return leftNames + rootName + rightNames
    else:
        return rightNames + rootName + leftNames

def getNode(tree, name):
    for n in tree.find_clades():
        if n.name == name:
            return n
    return None

def inorder_relative(tree, name):
    newtree = copy.deepcopy(tree)
    #node = newtree.find_any(name)
    node = getNode(newtree, name)
    if node:
        try:
            newtree.root_with_outgroup(node)
        except TypeError:
            print("Node %s" %node.name)
            print("Ancestor %s" % getParent(newtree, node).name)
            sys.exit(1)
        return inorder(newtree.root)
    else:
        sys.stderr.write("Cannot find node *%s* in tree using find_any\n" %name)
        return inorder(tree.root) 

def alignInternalNodes(tree):
    newtree = copy.deepcopy(tree)
    for clade in newtree.get_nonterminals():
        #continue #HACK
        #if clade.name != 'reference':#HACK
        #    continue#HACK
        children = clade.clades
        if len(children) != 2:
            sys.stderr.write("Tree %s is not a binary tree format! Binary tree is required. Please check.\n" %clade.name)
        selfleaf = Clade(branch_length=0, name=clade.name)
        newchildren=[children[0], selfleaf, children[1]]
        clade.clades = newchildren
    return newtree

def getLeaves(tree):
    return [ leaf.name for leaf in tree.get_terminals() ]

def getNodes(tree):
    return [ n.name for n in tree.find_clades() ]

def iterAllClades(root):
    allClades = []
    names = getNodes(root)
    children = []
    for name in names:
        if name:
            children.append(name)
    allClades.append(children)
    for clade in root.clades:
        allClades.extend( iterAllClades(clade) )
    return allClades

#=== get neighboring nodes ===
def getParent(tree, clade):
    for c in tree.find_clades():
        if clade in c.clades:
            return c
    return None

def getNeighbors(tree, name):
    #Return: 1 parent (if not root), 3 closest leaves, 1 leave that is a little further away
    neighbors = []
    allSortedNames = inorder_relative(tree, name)
    leaves =[c.name for c in tree.get_terminals()]
    sortedNames = []#leaves only
    for n in allSortedNames:
        if n in leaves:
            sortedNames.append(n)

    node = tree.find_any(name)
    parent = getParent(tree, node)
    if parent:
        neighbors.append(parent.name)
    
    currindex = 0
    for i, n in enumerate(sortedNames):
        if len(neighbors) >= 5:
            currindex = i
            break
        if n != name or (parent and n != parent.name):
            neighbors.append(n)
    if currindex < len(sortedNames):
        middleindex = (len(sortedNames) + currindex)/2
        neighbors.append(sortedNames[int(middleindex)])
    return neighbors


###### ETE2 #####
def my_layout(node):
    from ete2 import AttrFace, faces
    if node.is_leaf():
        faces.add_face_to_node(AttrFace("name"), node, column=0, position="aligned")
    else:
        node.img_style["size"] = 0
        node.img_style["fgcolor"] = "#000000"
        #node.img_style["shape"] = "circle"

def drawTree(nwfile, outfile):
    from ete2 import Tree, TreeStyle, TextFace
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    ts.branch_vertical_margin = 12.75
    ts.orientation = 1
    titleFace = TextFace("Phylogenetic Tree", fsize=18, fgcolor="white")
    titleFace.margin_top = 15
    ts.title.add_face(titleFace, column=1)

    t = Tree(nwfile)
    t.render(outfile, tree_style=ts)
    #t.render(outfile, w=183, units="mm", tree_style=ts)

#============ HAL RELATED ==========
def checkHalTree(halfile, outdir, options):
    treefile = os.path.join(outdir, "haltree.nw")
    system("halStats --tree %s > %s" %(halfile, treefile))
    tree = Phylo.read(treefile, "newick")
    options.treeFile = treefile
    options.tree = tree

#def getOrderFromTree(options):
#    if options.tree:
#        orders = inorder(options.tree, reverse=True)
#        if not options.genomes:
#            options.genomes = orders

def getProperName_tree(tree, properName):
    from hal.assemblyHub.assemblyHubCommon import getProperName
    for node in tree.find_clades():
        if node.name:
            node.name = getProperName(node.name, properName)
    return

def drawTreeWtInternalNodesAligned(tree, outdir, properName):
    modifiedTree = alignInternalNodes(tree)  
    leaves = getLeaves(modifiedTree)
    getProperName_tree(modifiedTree, properName)
    modifiedTreeFile = os.path.join(outdir, "hubTreeModified.nw")
    Phylo.write(modifiedTree, modifiedTreeFile, 'newick')
    treeFig = os.path.join(outdir, "hubTree.png")
    try:
        #drawTree(modifiedTreeFile, treeFig) #HACK
        return treeFig, leaves
    except RuntimeError:
        sys.stderr.write("Cannot draw the tree, so the tree will be missing on the config page. Perhaps ETE2 was not installed?\n")
        return None, None

