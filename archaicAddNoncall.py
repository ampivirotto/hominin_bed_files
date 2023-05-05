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
import sys

def REF(chrnum):
    """
    get reference genome
    """
    refDict = {}
    genome = pysam.FastaFile('/mnt/z/hg19_reference/chr{}.fa.gz'.format(chrnum))
    #import pdb; pdb.set_trace()
    for x in range(0, genome.lengths[0]):
        pos = genome.fetch(genome.references[0], start=x, end=x+1)
        refDict[x] = pos
    return refDict


def main(chrnum, vcffname):
    ## create new vcf output file
    outfname = vcffname.rstrip('.vcf.gz') + '_nonCallRegionAdded.vcf'
    New = open(outfname,"w")

    ref = REF(chrnum)
    print(len(ref))
    exit()

    with gzip.open(vcffname) as f:
        prevPos = 0
        for line in f:
            line = line.decode('ASCII')

            if line.startswith('#'):
                New.write(line)
            else:
                lineSplit = line.split('\t')
                curPos = lineSplit[1]



if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
