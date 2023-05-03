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

def VCFfilter():
    ##  bcftools view -R (coding region file) (original VCF file of Gorilla/Archaic) -i 'TYPE="snp"' | bcftools sort -o (name of new output VCF) -Oz

def index():
    ## index using bcftools

def zipUP():
    ## bgzip using bcftools

def main(location):
    ## split the large single wgs file into chromosomes (gorilla and pongo)

    ## iterate through every chromosome
    for file in os.listdir(location):
        if file.endswith('.vcf.gz'):
            ## filter out snps and sort the files
            VCFfilter()

            if not archaic:
                addNonCall()
            else:
                addArchaicMissing()


if __name__ == '__main__':
    ## folder name
    main(sys.argv[1])
