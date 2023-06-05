import pandas as pd
import numpy as np
import gzip
import sys

missing_stat  = float(sys.argv[1]) #system input for missingness
chrome = sys.argv[2]
names = sys.argv[3]

path = 'allSpeciesVCFs/chr{}_cds_primates.vcf.gz'.format(chrome)
vcf = gzip.open(path,'r')
idfile = open(names,'r')

def makeIDDictionary(line,file): #makes a dictionary of all the ids and the species they are
    indDict = {}
    specDict = {}

    df = pd.read_csv(file, sep = '\t', header = None)

    lineSplit = line.split('\t')

    for i in range(9, len(lineSplit)):
        idName = lineSplit[i]
        row = df[df[0] == idName]
        if len(row)>0:
            species = list(row[1])[0]

        indDict[i] = species

        if species not in specDict:
            specDict[species] = [i]
        else:
            specDict[species].append(i)

    return indDict, specDict

for line in vcf:   #getting the ids add putting in a list
    line = line.decode('ASCII')
    if line[:3] == '#CH':
        no,Go = makeIDDictionary(line,idfile)
        break


def makeCSV(vcf,ids,missing,chrom):
    output = open('AncesterStates{}.csv'.format(chrom))
    output.write('Chromosome,Position,MH Nuke,Primate Nuke'+'\n')
    for line in vcf:
        line = line.decode('ASCII')

        if line[0] != '#':
            MHFreq = {}
            PrimateFreq = {}
            lst = line.split() #list of values in line
            chrome = lst[0] #chrom in line
            pos = lst[1] #position in line
            nuke_ref = lst[3] #ref in line
            nuke_alt = lst[4] #alt in line
            Tot = 0
            ## figure out how many alt alleles there are
            triallelic = False
            quadallelic = False
            if len(nuke_alt.split(',')) == 2:
                triallelic = True
                nuke_alt, altBP2 = nuke_alt.split(',')
            elif len(nuke_alt.split(',')) == 3:
                quadallelic = True
                nuke_alt, altBP2, altBP3 = nuke_alt.split(',')


            for species in ids: #go species by species

                use = ids[species] #get list of ids
                for index in use: #go id by id for each species
                    data = lst[index][:3] #data for eact id
                    if species in ['Pan','Pongo','Gorilla']: #puts in data depending on species (just ancestors)
                        if len(data) < 2:
                            PrimateFreq['NA'] = PrimateFreq.get('NA',0) + 2
                        else:
                            PrimateFreq[nuke_ref] = PrimateFreq.get(nuke_ref,0) + data.count('0')
                            PrimateFreq[nuke_alt] = PrimateFreq.get(nuke_alt,0) + data.count('1')
                            PrimateFreq['NA'] = PrimateFreq.get('NA',0) + data.count('.')
                            if triallelic:
                                PrimateFreq[altBP2] = PrimateFreq.get(altBP2,0) + data.count('2')
                            if quadallelic:
                                PrimateFreq[altBP2] = PrimateFreq.get(altBP2,0) + data.count('2')
                                PrimateFreq[altBP3] = PrimateFreq.get(altBP3,0) + data.count('3')
                        Tot += 2
                    elif species == 'mHuman':
                        if len(data) < 2:
                            MHFreq['NA'] = MHFreq.get('NA',0) + 2
                        else:
                            MHFreq[nuke_ref] = MHFreq.get(nuke_ref,0) + data.count('0')
                            MHFreq[nuke_alt] = MHFreq.get(nuke_alt,0) + data.count('1')
                            MHFreq['NA'] = MHFreq.get('NA',0) + data.count('.')
                            if triallelic:
                                MHFreq[altBP2] = MHFreq.get(altBP2,0) + data.count('2') #data for all species
                            if quadallelic:
                                MHFreq[altBP2] = MHFreq.get(altBP2,0) + data.count('2')
                                MHFreq[altBP3] = MHFreq.get(altBP3,0) + data.count('3')
                        Tot += 2

            tot_missing = (MHFreq['NA'] + PrimateFreq['NA'])/Tot
            MHHighest = max(MHFreq)
            PrimateHighest = max(PrimateFreq)
            if tot_missing  < missing and MHHighest != 'NA':
                if MHHighest != PrimateHighest: 
                    output.write(','.join([chrome,pos,MHHighest,PrimateHighest])+'\n')
    output.close()



    

makeCSV(vcf,Go,missing_stat,chrome)


