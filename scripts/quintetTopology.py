'''


quintetTopology.py version 1.0 copyright (C) Joel Sharbrough 2020

USAGE:

	python quintetTopology.py trees.fofn outgroupID paternalDiploidID maternalDiploidID polyploidID > summaryFile.txt

Notes:

	-trees.fofn is a file containing file names for all the trees for which you want to test quintet topology on
	-Species IDs need to be unique identifiers by which genes from that species can be positively identified


'''

import dendropy
import sys


def quintetTopology(treeFile,outgroupID,paternalDiploidID,maternalDiploidID, polyploidID):
    tree = dendropy.Tree.get(path=treeFile,schema='newick',rooting='force-unrooted')
    fileSplit = treeFile.split('.')
    og = fileSplit[0]
    ogNum = fileSplit[1]
    totalLeaves = []
    for taxon in tree.taxon_namespace:
        totalLeaves.append(taxon.label)
        if outgroupID in taxon.label:
            outgroupLeaf = taxon.label
        elif maternalDiploidID in taxon.label:
            maternalDiploidLeaf = taxon.label
        elif paternalDiploidID in taxon.label:
            paternalDiploidLeaf = taxon.label
    tree.prune_taxa_with_labels([outgroupLeaf],update_bipartitions=True,suppress_unifurcations=False)
    treeFilter = False
    failingReasons = []
    rootPlacement = False
    for node in tree:
        subTree = node.leaf_nodes()
        currLeaves = []
        for leaf in subTree:
            currLeaves.append(leaf.taxon.label)
        if len(currLeaves) == 2:
            currSpeciesLeafDict,numSpecies = checkLeaves(currLeaves,[outgroupID,paternalDiploidID,maternalDiploidID,polyploidID])
            if currSpeciesLeafDict[maternalDiploidID] > 0 and currSpeciesLeafDict[paternalDiploidID] > 0:
                treeFilter = True
                failingReasons.append('Diploids group')
            elif currSpeciesLeafDict[polyploidID] > 1:
                treeFilter = True
                failingReasons.append('Tetraploids group')
    if treeFilter == True:
        sys.stdout.write(og + '\t' + ogNum + '\tFAILED\tFAILED')
        for item in failingReasons:
            sys.stdout.write('\t' + item)
        sys.stdout.write('\n')
    else:
        tree = dendropy.Tree.get(path=treeFile,schema='newick',rooting='force-unrooted')
        outgroupNode = tree.find_node_with_taxon_label(outgroupLeaf)
        tree.to_outgroup_position(outgroupNode, update_bipartitions=False)
        for node in tree:
            subTree = node.leaf_nodes()
            currLeaves = []
            for leaf in subTree:
                currLeaves.append(leaf.taxon.label)
            if len(currLeaves) == 3:
                treeFilter = True
                if maternalDiploidLeaf not in currLeaves:
                    rootPlacement = maternalDiploidID
                elif paternalDiploidLeaf not in currLeaves:
                    rootPlacement = paternalDiploidID
                else:
                    rootPlacement = polyploidID
    if treeFilter == True and rootPlacement != False:
        sys.stdout.write(og + '\t' + ogNum + '\tPASSED\tFAILED\t' + rootPlacement + '\n')
    elif treeFilter == False:
        sys.stdout.write(og + '\t' + ogNum + '\tPASSED\tPASSED\n')
        
    

def checkLeaves(leaves,speciesIDList):
    speciesLeafDict = {}
    for speciesID in speciesIDList:
        speciesLeafDict[speciesID] = 0
    for leaf in leaves:
        currSpecies = False
        for speciesID in speciesIDList:
            if speciesID in leaf:
                currSpecies = speciesID
                break
        if currSpecies != False:
            speciesLeafDict[currSpecies] += 1
    numSpecies = 0
    for species in speciesIDList:
        if speciesLeafDict[species] > 0:
            numSpecies += 1
    return speciesLeafDict,numSpecies

def midpointRoot(treeFile):
    tree = dendropy.Tree.get(path = treeFile,schema='newick')
    tree.reroot_at_midpoint(update_bipartitions=False)
    rootedTree = tree.as_string(schema='newick')
    return rootedTree

def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict


helpStatement='\n\n\nquintetTopology.py version 1.0 copyright (C) Joel Sharbrough 2020\n\nUSAGE:\n\n\tpython quintetTopology.py trees.fofn outgroupID paternalDiploidID maternalDiploidID polyploidID > summaryFile.txt\n\nNotes:\n\n\t-trees.fofn is a file containing file names for all the trees for which you want to test quintet topology on\n\t-Species IDs need to be unique identifiers by which genes from that species can be positively identified\n\n\n'
if len(sys.argv) == 1:
	sys.stderr.write(helpStatement)
elif sys.argv[1] == 'help' or sys.argv[1] == 'h':
	sys.stderr.write(helpStatement)
elif len(sys.argv) != 6:
	sys.stderr.write(helpStatement)
else:
	infile = open(sys.argv[1])
	sys.stdout.write('OGID\tsubTreeNumber\tTopology?\tRooting?\tReasons\n')
	for line in infile:
    	realLine = line
	    while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
    	    realLine = realLine[0:-1]
	    quintetTopology(realLine,sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) #treeFile,outgroupID, paternalDiploid, maternalDiploid, polyploid