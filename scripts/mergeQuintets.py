'''


mergeQuintets.py version 1.0 copyright (C) Joel Sharbrough 2020

USAGE:

	python mergeQuintets.py phylologs.txt syntelogs.txt outgroupID paternalDiploidID paternalTetraploidID maternalTetraploidID maternalDiploidID > summaryFile.txt

Notes:
	
	-phylologs.txt is a tab-delimited text file containing the gene names for each single copy quintet produced by phylogenetics
	-syntelogs.txt is a tab-delimited text file containing the gene names for each single copy quintet produced by synteny analysis
	-Species IDs need to be unique identifiers by which genes from that species can be positively identified


'''

import sys

def mergeQuintets(phyloLogs,synteLogs,outgroupID,paternalDiploidID,paternalTetraploidID,maternalTetraploidID,maternalDiploidID):
    quintetSplit = phyloLogs.split('_')
    infile = open(phyloLogs,'r')
    phyloDict = {}
    phyloGeneList = []
    phyloGroupList = []
    for line in infile:
        if line[0] != '#':
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split('\t')
            phyloGeneList += lineSplit[1:]
            phyloGroupList.append(lineSplit[0])
            phyloDict[lineSplit[0]] = lineSplit[1:]
    infile.close()
    sys.stderr.write('phyloList contains ' + str(len(phyloGroupList)) + ' singleton genes\n')
    infile = open(synteLogs,'r')
    synteGeneList = []
    synteGroupList = []
    synteDict = {}
    for line in infile:
        if line[0] != '#':
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split('\t')
            synteGeneList += lineSplit[1:]
            synteGroupList.append(lineSplit[0])
            synteDict[lineSplit[0]] = lineSplit[1:]
    sys.stderr.write('synteList contains ' + str(len(synteGroupList)) + ' singleton genes\n')
    infile.close()
    i = 0
    outfile1 = open(quintetSplit[0] + '.identicalQuintets.txt','w')
    outfile1.write('#OGID\toutgroup\tpaternalDiploid\tpaternalTetraploid\tmaternalTetraploid\tmaternalDiploid\n')
    outfile2 = open(quintetSplit[0] + '.overlappingQuintets.txt','w')
    outfile2.write('#OGID\toutgroup\tpaternalDiploid\tpaternalTetraploid\tmaternalTetraploid\tmaternalDiploid\n')
    outfile3 = open(quintetSplit[0] + '.phylogeneticSpecificQuintets.txt','w')
    outfile3.write('#OGID\toutgroup\tpaternalDiploid\tpaternalTetraploid\tmaternalTetraploid\tmaternalDiploid\n')
    outfile4 = open(quintetSplit[0] + '.syntenySpecificQuintets.txt','w')
    outfile4.write('#OGID\toutgroup\tpaternalDiploid\tpaternalTetraploid\tmaternalTetraploid\tmaternalDiploid\n')
    sys.stdout.write('#OGID\tOG Merge Status\tphyloGroup(s)\tsynteGroup(s)\n')
    phyloSpecific = []
    synteSpecific = []
    phyloMerged = []
    synteMerged = []
    phyloOverlap = []
    synteOverlap = []
    numOverlaps = 0
    for phyloGroup in phyloGroupList:
        match = False
        j = 0
        currPhyloGenes = phyloDict[phyloGroup]
        '''if len(currPhyloGenes) == 5:  #CHANGE THIS FOR COTTTON
            currPhyloGenes = currPhyloGenes[0:-1]'''
        while match == False and j < len(synteGroupList):
            currSynteGroup = synteGroupList[j]
            if currPhyloGenes == synteDict[currSynteGroup]:
                match = currSynteGroup
            j += 1
        if match != False:
            phyloMerged.append(phyloGroup)
            synteMerged.append(match)
            outfile1.write('OG' + str(i) + '\t' + phyloGroup + '\t' + match + '\t')
            sys.stdout.write('OG' + str(i) + '\tIdentical\t' + phyloGroup + '\t' + match + '\n')
            for seq in phyloDict[phyloGroup]:
                outfile1.write(seq + '\t')
            outfile1.write('\n')
        else:
            overlap = False
            for seq in phyloDict[phyloGroup]:
                if seq in synteGeneList:
                    overlap = True
            if overlap == False:
                outfile3.write('OG' + str(i) + '\t' + phyloGroup + '\tN/A\t')
                sys.stdout.write('OG' + str(i) + '\tPhylo-specific\t' + phyloGroup + '\tN/A\n')
                for seq in phyloDict[phyloGroup]:
                    outfile3.write(seq + '\t')
                outfile3.write('\n')
                phyloSpecific.append(phyloGroup)
            else:
                numOverlaps += 1
                phyloGenes = phyloDict[phyloGroup]
                synteGenes = []
                phyloGroups = [phyloGroup]
                synteGroups = []
                k = 0
                while k < len(phyloGenes):
                    for synteGroup in synteGroupList:
                        if phyloGenes[k] in synteDict[synteGroup] and synteGroup not in synteGroups:
                            synteGroups.append(synteGroup)
                            for gene in synteDict[synteGroup]:
                                if gene not in synteGenes:
                                    synteGenes.append(gene)
                    l = 0
                    while l < len(synteGenes):
                        for phyloGroup2 in phyloGroupList:
                            if synteGenes[l] in phyloDict[phyloGroup2] and phyloGroup2 not in phyloGroups:
                                phyloGroups.append(phyloGroup2)
                                for gene in phyloDict[phyloGroup2]:
                                    if gene not in phyloGenes:
                                        phyloGenes.append(gene)
                        l += 1
                    k += 1
                outfile2.write('OG' + str(i) + '\t' + phyloGroups[0])
                sys.stdout.write('OG' + str(i) + '\tConflicting\t' + phyloGroups[0])
                if len(phyloGroups) > 1:
                    for pg in phyloGroups[1:]:
                        outfile2.write(',' + pg)
                        sys.stdout.write(',' + pg)
                outfile2.write('\t' + synteGroups[0])
                sys.stdout.write('\t' + synteGroups[0])
                if len(synteGroups) > 1:
                    for sg in synteGroups[1:]:
                        outfile2.write(',' + sg)
                        sys.stdout.write(',' + sg)
                outfile2.write('\t')
                sys.stdout.write('\n')
                phyloOverlap += phyloGroups
                synteOverlap += synteGroups
                outgroupSeqs = []
                paternalDiploidSeqs = []
                paternalTetraploidSeqs = []
                maternalDiploidSeqs = []
                maternalTetraploidSeqs = []
                allGenes = phyloGenes + synteGenes
                for seq in allGenes:
                    if outgroupID in seq and seq not in outgroupSeqs:
                        outgroupSeqs.append(seq)
                    elif paternalDiploidID in seq and seq not in paternalDiploidSeqs:
                        paternalDiploidSeqs.append(seq)
                    elif paternalTetraploidID in seq and seq not in paternalTetraploidSeqs:
                        paternalTetraploidSeqs.append(seq)
                    elif maternalTetraploidID in seq and seq not in maternalTetraploidSeqs:
                        maternalTetraploidSeqs.append(seq)
                    elif maternalDiploidID in seq and seq not in maternalDiploidSeqs:
                        maternalDiploidSeqs.append(seq)
                if len(outgroupSeqs) == 1:
                    outfile2.write(outgroupSeqs[0])
                elif len(outgroupSeqs) > 1:
                    outfile2.write(outgroupSeqs[0])
                    for outSeq in outgroupSeqs[1:]:
                        outfile2.write(',' + outSeq)
                outfile2.write('\t')
                if len(paternalDiploidSeqs) == 1:
                    outfile2.write(paternalDiploidSeqs[0])
                elif len(paternalDiploidSeqs) > 1:
                    outfile2.write(paternalDiploidSeqs[0])
                    for paternalDiploidSeq in paternalDiploidSeqs[1:]:
                        outfile2.write(',' + paternalDiploidSeq)
                outfile2.write('\t')
                if len(paternalTetraploidSeqs) == 1:
                    outfile2.write(paternalTetraploidSeqs[0])
                elif len(paternalTetraploidSeqs) > 1:
                    outfile2.write(paternalTetraploidSeqs[0])
                    for paternalTetraploidSeq in paternalTetraploidSeqs[1:]:
                        outfile2.write(',' + paternalTetraploidSeq)
                outfile2.write('\t')        
                if len(maternalTetraploidSeqs) == 1:
                    outfile2.write(maternalTetraploidSeqs[0])
                elif len(maternalTetraploidSeqs) > 1:
                    outfile2.write(maternalTetraploidSeqs[0])
                    for maternalTetraploidSeq in maternalTetraploidSeqs[1:]:
                        outfile2.write(',' + maternalTetraploidSeq)
                outfile2.write('\t')
                if len(maternalDiploidSeqs) == 1:
                    outfile2.write(maternalDiploidSeqs[0])
                elif len(maternalDiploidSeqs) > 1:
                    outfile2.write(maternalDiploidSeqs[0])
                    for maternalDiploidSeq in maternalDiploidSeqs[1:]:
                        outfile2.write(',' + maternalDiploidSeq)
                outfile2.write('\n')
        i += 1
    for synteGroup in synteGroupList:
        overlap = False
        for seq in synteDict[synteGroup]:
            for pg in phyloGroupList:
                if seq in phyloDict[pg]:
                    overlap = True
        if overlap == False:
            synteGenes = synteDict[synteGroup]
            outfile4.write('OG' + str(i) + '\tN/A\t' + synteGroup + '\t')
            sys.stdout.write('OG' + str(i) + '\tSynteny-specific\tN/A\t' + synteGroup + '\n')
            synteSpecific.append(synteGroup)
            for seq in synteDict[synteGroup]:
                outfile4.write(seq + '\t')
            outfile4.write('\n')
            i += 1
    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()
    sys.stderr.write('In total there were:\n\t' + str(len(phyloMerged)) + ' identical quintets\n\t' + str(len(phyloOverlap)) + ' phylogenetic groups and ' + str(len(synteOverlap)) + ' syntenic groups that together formed ' + str(numOverlaps) + ' non-identical, overlapping groups\n\t' + str(len(phyloSpecific)) + ' phylo-specific quintets\n\t' + str(len(synteSpecific)) + ' synteny-specific quintets\n')
            
helpStatement='\n\n\nmergeQuintets.py version 1.0 copyright (C) Joel Sharbrough 2020\n\nUSAGE:\n\n\tpython mergeQuintets.py phylologs.txt syntelogs.txt outgroupID paternalDiploidID paternalTetraploidID maternalTetraploidID maternalDiploidID > summaryFile.txt\n\nNotes:\n\n\t-Each species ID is a unique identifier by which genes from that species can be positively identified\n\n\n'
if len(sys.argv) == 1:
	sys.stderr.write(helpStatement)
elif sys.argv[1] == 'help' or sys.argv[1] == 'h':
	sys.stderr.write(helpStatement)
elif len(sys.argv) != 8:
	sys.stderr.write(helpStatement)
else:
	mergeQuintets(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])