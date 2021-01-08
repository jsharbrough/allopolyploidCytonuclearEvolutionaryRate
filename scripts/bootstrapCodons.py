import sys
import random

def bootstrapCodons(pamlAlignment):
    infile = open(pamlAlignment,'r')
    lines = infile.readlines()
    infile.close()
    totalLine = lines[0]
    while totalLine[-1] == '\t' or totalLine[-1] == '\n' or totalLine[-1] == '\r':
        totalLine = totalLine[0:-1]
    totalLineSplit = totalLine.split(' ')
    totalLength = int(totalLineSplit[-1])
    codonDict = {}
    i = 2
    speciesList = []
    while i < len(lines):
        j = 0
        currSeq = lines[i]
        while currSeq[-1] == '\t' or currSeq[-1] == '\r' or currSeq[-1] == '\n':
            currSeq = currSeq[0:-1]
        currSpecies = lines[i-1]
        while currSpecies[-1] == '\t' or currSpecies[-1] == '\r' or currSpecies[-1] == '\n':
            currSpecies = currSpecies[0:-1]
        codonList = []
        while j < totalLength:
            codonList.append(currSeq[j] + currSeq[j+1] + currSeq[j+2])
            j += 3
        codonDict[currSpecies] = codonList
        speciesList.append(currSpecies)
        i += 2
    bootCodons = []
    for num in range(totalLength/3):
        bootCodons.append(random.choice(range(totalLength/3)))
    sys.stdout.write(' 5 ' + str(totalLength) + '\n')
    for species in speciesList:
        sys.stdout.write(species + '\n')
        currCodons = codonDict[species]
        bootSeq = ''
        for codonNum in bootCodons:
            bootSeq += currCodons[codonNum]
        sys.stdout.write(bootSeq +  '\n')
bootstrapCodons(sys.argv[1])