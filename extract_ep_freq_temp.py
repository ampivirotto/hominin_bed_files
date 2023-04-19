#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     11/04/2023
# Copyright:   (c) tuk32868 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import gzip
import pandas as pd
import numpy as np


def readVCF(vcffn, indLoc):
    with open(vcffn) as f:
        for line in f:
            if not line.startswith('#'):
                readLine(line, indDict, spD)
            elif line.startswith('#C'):
                indDict, spD = makeIndDictionary(line, indLoc)


def readGZIPvcf(vcffn, indLoc):
    with gzip.open(vcffn) as f:
        for line in f:
            line = line.decode('ASCII')
            if not line.startswith("#"):
                readLine(line, indDict, spD)
            elif line.startswith('#C'):
                indDict, spD = makeIndDictionary(line, indLoc)


def makeIndDictionary(line, location):
    indDict = {}
    specDict = {}

    df = pd.read_csv('ep/primate_branch_ids.txt', sep = '\t', header = None)

    lineSplit = line.split('\t')

    for i in range(9, len(lineSplit)):
        idName = lineSplit[i]
        species = list(df[df[0] == idName][1])[0]

        indDict[i] = species

        if species not in specDict:
            specDict[species] = [i]
        else:
            specDict[species].append(i)

    return indDict, specDict


def checkFreq(ids, line, secondAlt = False):
    missing, ref, alt = 0, 0, 0
    if secondAlt:
        alt2 = 0

    for idi in ids:
        gt = line[idi][:3]

        if len(gt) == 1:
            missing += 2
        else:
            alt += gt.count('1')
            ref += gt.count('0')
            missing = gt.count('.')
            if secondAlt:
                alt2 = gt.count('2')

    if secondAlt:
        return ref/len(ids), alt/len(ids), missing/len(ids), alt2/len(ids)
    return ref/len(ids), alt/len(ids), missing/len(ids)

def readLine(line, indDict):
    splitLine = line.strip('\n').split('\t')

    pos = splitLine[1]
    chrnum = splitLine[0]
    ref = splitLine[3]
    alt = splitLine[4]

    mhREF, mhALT, mhMissing = checkFreq(specs['mHuman'], line)

    if mhMissing > 0.5:
        return





def main(vcffn):

    if vcffn.endswith('.gz'):
        readGZIPvcf(vcffn)
    else:
        readVCF(vcffn)

if __name__ == '__main__':
    main()
