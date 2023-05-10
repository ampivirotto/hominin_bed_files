

import pandas as pd
import numpy as np
import gzip
import pickle
import sys

def turnBedintoDict(bedfn, log):
    bedDict = {}
    with open(bedfn) as f:
        for line in f:
            row = line.split('\t')
            for pos in np.arange(int(row[1]), int(row[2])):
                if 'chr' in row[0]:
                    chrnum = row[0].lstrip('chr')
                else:
                    chrnum = row[0]
                if pos in bedDict:
                    log.write('Duplicate: ' + str(pos) + '\n')
                bedDict[int(pos)] = [str(chrnum), str(pos)]
    return bedDict

def addREF(mainDict, refdict, colCount):
    colCount += 1
    for item in mainDict.keys():
        ref = refdict[item]

        mainDict[item].append(ref)

        assert len(mainDict[item]) == colCount
    return mainDict, colCount

def addVCFcols(mainDict, columnCount):
    columnCount += 1
    for item in mainDict.keys():
        mainDict[item].append('.')
        assert len(mainDict[item]) == columnCount
    return mainDict, columnCount

def addVCFindividuals(mainDict, numInds, colCount):
    colCount += numInds
    for item in mainDict.keys():
        for indN in range(numInds):
            mainDict[item].append('./.')
        assert len(mainDict[item]) == colCount
    return mainDict, colCount

def main(chromNum, ncBED, vcffname):
    ## create new vcf output file
    outfname = vcffname.rstrip('.vcf.gz') + '_nonCallRegionAdded.vcf'
    New = open(outfname,"w")

    ## create log file
    log = open(vcffname.rstrip('.vcf.gz') + "_log.txt", "w")

    ## read in bed file of noncallable regions
    NonCall = turnBedintoDict(ncBED, log)
    colCount = 2

    ## read in dictionary of ancestor
    dic = pd.read_pickle(r'./hg19_ancestor_pickle/grch37_poskey_baseval_dic_{}_based.p'.format(chromNum))

    ## add id col
    workingDict, colCount = addVCFcols(NonCall, colCount)

    ## select rows for noncallable
    workingDict, colCount = addREF(NonCall, dic, colCount)

    ## turn dfs into a row df for each site in vcf format
    for x in range(5):
        workingDict, colCount = addVCFcols(workingDict, colCount)


    missingSites = len(workingDict)
    log.write('Total Noncallable Sites in CDS region: ' + str(missingSites) + '\n')
    #print(missingSites)

    vcfSites = 0
    foundNC = 0
    with gzip.open(vcffname, 'r') as file:
        for line in file:
            line = line.decode('ASCII')
            if line[0] != '#':
                temp = line.strip('\n').split("\t")
                if int(temp[1]) in workingDict:
                    log.write(line)
                    foundNC += 1
                else:
                    workingDict[int(temp[1])] = temp
                    vcfSites+=1
            elif line.startswith('#CHROM'):
                numInds = len(line.split('\t')) - 9
                workingDict, colCount = addVCFindividuals(workingDict, numInds, colCount)
                New.write(line)
            else:
                New.write(line)

    #print(len(workingDict))
    #print(missingSites + vcfSites)
    assert len(workingDict) == missingSites + vcfSites
    #print(vcfSites, foundNC)
    log.write('Total Number of Variant sites (in callable CDS): ' + str(vcfSites) +'\n')
    log.write('Total Number of Variant Sites (in noncallable CDS): ' + str(foundNC) + '\n')

    # --- Alexander was here ---
    for row in sorted(workingDict.keys()):
        new_line = '\t'.join(workingDict[row])
        New.write(new_line + '\n')

    New.close()

    return outfname

if __name__ == '__main__':
    ## python NAME.py CHROM BEDFILE VCFFILE
    outfn = main(sys.argv[1], sys.argv[2], sys.argv[3])
