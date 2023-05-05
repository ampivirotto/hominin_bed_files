import os
import subprocess
import re
import sys
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
def getChrome(file):
    to_search = r'chr\d+'
    all = re.search(to_search,file)
    if all == None:
        return 'all'
    return all[3:]

def VCFfilter(folder,file,chrome):
    if chrome == 'all':
        coding_region_file = './bedfiles/EP_CDS_no_CHR.merged.bed'
    else:
        coding_region_file = './bedfiles/'+ chrome +'.merged.bed'
    title = file.split('.')
    go = title.index('.vcf')
    title[go-1] = title[go-1] + '_filtered'
    output =  folder + '/' + '.'.join(title)
    command  = 'bcf view -R '+ coding_region_file +' ' + file + ''' -i 'TYPE = "snp"' | bcftools sort -o ''' + output + ' -Oz'
    ##  bcftools view -R (coding region file) (original VCF file of Gorilla/Archaic) -i 'TYPE="snp"' | bcftools sort -o (name of new output VCF) -Oz
    subprocess.run(command)
    return output

def index(file):
    ## index using bcftools
    command = 'bcftools index ' + file
    subprocess.run(command)

def zipUP(file):
    ## bgzip using bcftools
    command = 'bgzip ' + file
    subprocess.run(command)

def main(location):
    ## split the large single wgs file into chromosomes (gorilla and pongo)

    ## iterate through every chromosome
    for file in os.listdir(location):
        if file.endswith('.vcf.gz'):
            ## filter out snps and sort the files
            chromosome = getChrome(file)
            done = VCFfilter(location,file,chromosome)
            index(done)
            zipUP(done)
            if not archaic:
                addNonCall()
            else:
                addArchaicMissing()


if __name__ == '__main__':
    ## folder name
    main(sys.argv[1])
