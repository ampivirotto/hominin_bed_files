
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     03/05/2023
# Copyright:   (c) tuk32868 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
import subprocess
import sys
import gzip
import shlex
import AddNoncall
import re
import ArchaicaddNonCall

def VCFfilter(folder,file,chrome):
    if chrome == 'all':
        coding_region_file = './bedfiles/EP_CDS_no_CHR.merged.bed'
    else:
        coding_region_file = './bedfiles/'+ chrome +'.merged.bed'
    newfile = file.rstrip('.vcf.gz') + '_filter.vcf.gz'
    output =  folder + '/' + newfile
    command  = 'bcf view -R '+ coding_region_file +' ' + file + ''' -i 'TYPE = "snp"' | bcftools sort -o ''' + output + ' -Oz'
    ##  bcftools view -R (coding region file) (original VCF file of Gorilla/Archaic) -i 'TYPE="snp"' | bcftools sort -o (name of new output VCF) -Oz
    cmd = shlex.split(command)
    subprocess.run(cmd)
    return output

def index(file):
    ## index using bcftools
    cmdLine = shlex.split('bcftools index {}'.format(file))
    subprocess.run(cmdLine)

def getChrome(file):
    to_search = r'chr\d+'
    all = re.search(to_search,file)
    if all == None:
        return 'all'
    return all[3:]

def zipUP(file):
    ## bgzip using bcftools
    command = 'bgzip ' + file
    command = shlex.split(command)
    subprocess.run(command)

def ReName(location,file):
    output = file.strip('.vcf.gz') + '_renamed_chr.vcf.gz'
    command = 'bcftools annotate --rename-chrs Chrome/chrName.txt -Oz -o ' + output
    cmdLine = shlex.split(command)
    subprocess.run(cmdLine)
    return output

def peakVCF(file):
    """
    identify if chromosome column is chr# or # format
    """
    if file.endswith('.gz'):
        with gzip.open(file) as f:
            for line in f:
                line = line.decode('ASCII')
                if not line.startswith("#"):
                    chrnum = line.split('\t')[0]
                    if 'chr' in chrnum:
                        return ''
                    else:
                        return 'chr'


def nonCall(location,file,chrom,bedfile):
    """
    add uncallable regions using AddNoncall.py program which takes three arguments: chromosome number, bedfile, vcffile
    """
    AddNoncall.main(chrom, bedfile, location + file)

def splitByChromosome(location):
    """
    split all vcf files into single chromosome vcf files
    """
    files = os.listdir(location)
    returnFiles = []

    for file in files:
        if file.endswith('.vcf.gz'):
            if not os.path.isfile(location + file + '.csi') or os.path.isfile(location + file + '.tbi'):
                index(location + file)
            chrform = peakVCF(location + file) #might change to getChrome()
            splitCMD = 'sh split.sh {} {} {}'.format(location, file, chrform)
            cmdLine = shlex.split(splitCMD)
            returnFiles.append(file)
            subprocess.run(cmdLine)
    return returnFiles

def main(location):
    ## split the large single wgs file into chromosomes (gorilla and pongo)
    location = sys.argv[1]
    if len(os.listdir(location)) < 22:
        new_Loc = location + '/WholeGenome'
        Files = splitByChromosome(new_Loc)
        for file in Files:
            command = 'mv ' + file + '~/hominin/' + location
            cmdLine = shlex.split(command)
            subprocess.run(cmdLine)



    ## iterate through every chromosome
    for file in os.listdir(location):
        if file.endswith('.vcf.gz'):
            ## filter out snps and sort the files
            chromosome = getChrome(file)
            output = VCFfilter(location,file,chromosome)
            if not archaic:
                nonCall()
            else:
                addArchaicMissing() #currently in the making


if __name__ == '__main__':
    ## folder name
    main(sys.argv[1])
