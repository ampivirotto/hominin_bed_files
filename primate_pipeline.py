
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
import archaicAddNoncall

def VCFfilter(folder,file,chrome):
    coding_region_file = './bedfiles/'+ chrome +'.merged.bed'
    newfile = file.rstrip('.vcf.gz') + '_filter.vcf.gz'
    output =  newfile
    command  = 'bcftools view -R {} {} -i \'TYPE = "snp"\' | bcftools sort -o {} -Oz'.format(coding_region_file, file, output)
    ##  bcftools view -R (coding region file) (original VCF file of Gorilla/Archaic) -i 'TYPE="snp"' | bcftools sort -o (name of new output VCF) -Oz
    cmd = shlex.split(command)
    subprocess.run(cmd)

    ## move file
    subfolder = folder + 'preFilter/'
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    moveFile(file, subfolder)
    if os.path.isfile(file + '.csi'):
        moveFile(file + '.csi', subfolder)
    elif os.path.isfile(file + '.tbi'):
        moveFile(file + '.tbi', subfolder)

    return output

def VCFarchFilter(folder, file, chrome):
    """
    filter based on region file and sort (NO filter based on snp status - removes reference sites)
    """
    coding_region_file = './bedfiles/'+ chrome +'.merged.bed'
    newfile = file.rstrip('.vcf.gz') + '_filter.vcf.gz'
    output =  newfile
    command  = 'bcftools view -R {} {} | bcftools sort -o {} -Oz'.format(coding_region_file, file, output)
    ##  bcftools view -R (coding region file) (original VCF file of Gorilla/Archaic) -i 'TYPE="snp"' | bcftools sort -o (name of new output VCF) -Oz
    cmd = shlex.split(command)
    subprocess.run(cmd)

    ## move file
    subfolder = folder + 'preFilter/'
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    moveFile(file, subfolder)
    if os.path.isfile(file + '.csi'):
        moveFile(file + '.csi', subfolder)
    elif os.path.isfile(file + '.tbi'):
        moveFile(file + '.tbi', subfolder)

    return output

def index(file):
    ## index using bcftools
    cmdLine = shlex.split('bcftools index {}'.format(file))
    subprocess.run(cmdLine)

def getChrome(file):
    to_search = r'chr\d+'
    all = re.findall(to_search,file) #re.findall(to_search,gzip.open(file).read())
    if len(all) == 0:
        pass #what to put of now chrom
    return all[0][3:]

def zipUP(file):
    ## bgzip using bcftools
    command = 'bgzip {}'.format(file)
    command = shlex.split(command)
    subprocess.run(command)

def moveFile(file, newFolder):
    command = 'mv {} {}'.format(file, newFolder)
    cmdLine = shlex.split(command)
    subprocess.run(cmdLine)

def ReName(location,file):

    ## rename the file from chr# to # format
    output = location + file.rstrip('.vcf.gz') + '_renamedCHR.vcf.gz'
    command = 'bcftools annotate --rename-chrs chrName.txt -Oz -o {} {}'.format(output, location + file)
    cmdLine = shlex.split(command)
    subprocess.run(cmdLine)

    ## move old file to a subfolder
    subfolder = location + 'chrVCF/'
    if not os.path.exists(subfolder):
            os.makedirs(subfolder)
    moveFile(location + file, subfolder)
    if os.path.isfile(location + file + '.csi'):
        moveFile(location + file + '.csi', subfolder)
    elif os.path.isfile(location + file + '.tbi'):
        moveFile(location + file + '.tbi', subfolder)
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
                        return '', chrnum[3:]
                    else:
                        return 'chr', chrnum

def ArchNonCall(location, vcffile, chrom):
    """
    add uncallable regions using archaicAddNoncall.py program which takes two arguments: chromosome number and vcffile
    """
    try:
        outfname = archaicAddNoncall.main(chrom, vcffile)
    except:
        outfname = ''
        print('ERROR in adding noncallable sites {}'.format(vcffile))

    ## move old file to a subfolder
    subfolder = location + 'postFilter/'
    if not os.path.exists(subfolder):
            os.makedirs(subfolder)
    moveFile(vcffile, subfolder)
    if os.path.isfile(vcffile + '.csi'):
        moveFile(vcffile + '.csi', subfolder)
    elif os.path.isfile(vcffile + '.tbi'):
        moveFile(vcffile + '.tbi', subfolder)
    return outfname

def nonCall(location,vcffile,chrom,bedfileLocation):
    """
    add uncallable regions using AddNoncall.py program which takes three arguments: chromosome number, bedfile, vcffile
    """
    bedfiles = os.listdir(bedfileLocation)
    tempBed = bedfiles[0].split('_')[:-1]
    bedfileName = "_".join(tempBed) + "_{}.bed".format(str(chrom))

    bedfile = bedfileLocation + bedfileName
    if not os.path.isfile(bedfile):
        print('error BED file for noncallable region not found for chromosome {}'.format(chrom))
        exit()
    try:
        newFile = AddNoncall.main(chrom, bedfile, vcffile)
    except:
        newFile = ''
        print('ERROR in adding noncallable sites {}'.format(vcffile))

    ## move old file to a subfolder
    subfolder = location + 'postFilter/'
    if not os.path.exists(subfolder):
            os.makedirs(subfolder)
    moveFile(vcffile, subfolder)
    if os.path.isfile(vcffile + '.csi'):
        moveFile(vcffile + '.csi', subfolder)
    elif os.path.isfile(vcffile + '.tbi'):
        moveFile(vcffile + '.tbi', subfolder)
    return newFile

def splitByChromosome(location):
    """
    split all vcf files into single chromosome vcf files
    """
    files = os.listdir(location)

    for file in files:
        if file.endswith('.vcf.gz'):
            if not os.path.isfile(location + file + '.csi') or os.path.isfile(location + file + '.tbi'):
                index(location + file)
            chrform, num = peakVCF(location + file) #getChrome(location + file) #might change to peakVcf()
            splitCMD = 'sh split.sh {} {} {}'.format(location, file, chrform)
            cmdLine = shlex.split(splitCMD)
            subprocess.run(cmdLine)
    return files

def main(location, bedfileLoc):
    ## split the large single wgs file into chromosomes (gorilla and pongo)
    if len(os.listdir(location)) < 22:
        Files = splitByChromosome(location)
        newloc = location + 'WholeGenome/'
        if not os.path.exists(newloc):
            os.makedirs(newloc)
        for file in Files:
            if os.path.isfile(location + file):
                command = 'mv {} {}'.format(location + file, newloc)
                cmdLine = shlex.split(command)
                subprocess.run(cmdLine)


    ## iterate through every chromosome
    curFiles = os.listdir(location)
    for file in curFiles:
        if file.endswith('.vcf.gz'):
            chromosome, num = peakVCF(location + file)
            #print(chromosome, num)

            if not 'Archaic' in location:
                ## filter out snps and sort the files
                if chromosome == '':
                    New_chrome_file = ReName(location,file)
                    index(New_chrome_file)
                    output = VCFfilter(location,New_chrome_file,num)
                    index(output)
                else:
                    output = VCFfilter(location,file,num)
                    index(output)

                newFileName = nonCall(location,output,num,bedfileLoc)
            else:
                if chromosome == '':
                    New_chrome_file = ReName(location,file)
                    index(New_chrome_file)
                    output = VCFarchFilter(location,New_chrome_file,num)
                    index(output)
                else:
                    vcffile = location + file
                    if not (os.path.isfile(vcffile + '.csi') or os.path.isfile(vcffile + '.tbi')):
                        index(vcffile)
                    output = VCFarchFilter(location,vcffile,num)
                    index(output)
                newFileName = ArchNonCall(location, output, num)

            if not newFileName == '':
                    zipUP(newFileName)
                    index(newFileName+'.gz')
            print('Finish {}'.format(location + file))


if __name__ == '__main__':
    ## folder name
    main(sys.argv[1], sys.argv[2])
    #main('/mnt/z/hominin_adaptation_vcfs/test_Pipeline/Archaics/Altai/', '')
    #os.chdir(__file__.rstrip("primate_pipeline.py"))
    #main('../Pongo/', '/mnt/z/hominin_adaptation_vcfs/callable_uncallable/Pongo_uncallable/')
