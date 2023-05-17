#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     04/05/2023
# Copyright:   (c) tuk32868 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import gzip
import pysam
import numpy as np
import sys

def REF(chrnum):
    """
    get reference genome
    """
    refDict = {}
    genome = pysam.FastaFile('./hg19_reference/chr{}.fa.gz'.format(chrnum)) #reading in file using pysam
    with open('./bedfiles/{}.merged.bed'.format(chrnum)) as f:
        for line in f:
            splitLine = line.strip('\n').split('\t')
            for x in range(int(splitLine[1]), int(splitLine[2])+1):
                pos = genome.fetch(genome.references[0], start=x, end=x+1) #fetching the position
                refDict[x] = pos
    return refDict

def makeLine(chrnum, pos, ref, numInds):
    info = '.'
    if ref.islower():
        info = 'hg19REF=lowProbability;'
    cols = [chrnum, pos, '.', ref, '.', '.', '.', info, '.']  ## chr, pos, id, ref, alt, qual, filter, info, format
    assert len(cols) == 9

    for indn in range(numInds):
        cols.append('./.')

    return '\t'.join(cols) + '\n'

def main(chrnum, vcffname):
    ## create new vcf output file
    outfname = vcffname.rstrip('.vcf.gz') + '_nonCallRegionAdded.vcf'
    New = open(outfname,"w")

    ref = REF(chrnum)

    with gzip.open(vcffname) as f:
        prevPos = min(ref.keys())
        for line in f:
            try:
                line = line.decode('ASCII')

                if line.startswith('#CHROM'):
                    numInds = len(line.split('\t')) - 9
                    New.write(line)
                elif line.startswith('#'):
                    New.write(line)
                else:
                    lineSplit = line.split('\t')
                    curPos = int(lineSplit[1])
                    if curPos > prevPos:
                        for x in np.arange(prevPos, curPos):
                            if x in ref.keys():
                                if ref[x] != 'N':
                                    New.write(makeLine(chrnum, str(x), ref[x], numInds))
                    prevPos = curPos + 1
                    New.write(line)
            except:
                print('here')
    return outfname



if __name__ == '__main__':
    outfname = main(sys.argv[1], sys.argv[2])
    #outfname = main('22', 'Archaics/chr22_mq25_mapab100_filtered.vcf.gz')
