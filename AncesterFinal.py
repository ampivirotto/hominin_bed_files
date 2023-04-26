import pandas as pd
import numpy as np
import gzip
import sys

missing_stat  = float(sys.argv[1]) #system input for missingness
vcf = sys.argv[2]
names = sys.argv[3]

vcf = gzip.open(vcf,'r')
idfile = open(names,'r')

def makeIndDictionary(line,file):
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
        no,Go = makeIndDictionary(line,idfile)
        break


def makeTXT(vcf,ids,missing):
    cols = ['Chromosome','Position','Ref','Alt', 'Alt2', 'Alt3']
    for species in ids.keys():
        for x in ['ref', 'alt', 'alt2', 'alt3', 'missing']:
            cols.append(species + "_" + x)
    output_df_dic = {}
    #above is a dataframe just to easily organize the data
    state_df_dic = {} #Ensembl will be added later

    with open(vcf,'r') as vcf:
        for line in vcf:
            line = line.decode('ASCII')

            if line[0] != '#':
                add = [] #lest to add to freq file
                state_add = [] #list to add to state file
                lst = line.split() #list of values in line
                chrome = lst[0] #chrom in line
                pos = lst[1] #position in line
                nuke_ref = lst[3] #ref in line
                nuke_alt = lst[4] #alt in line
                ancest_missing, ancest_tot, ancest_alt, ancest_ref = 0, 0, 0, 0 #initiate for acestors
                an_alt2, an_alt3 = 0, 0
                state_add.extend([chrome,pos])

                ## figure out how many alt alleles there are
                triallelic = False
                quadallelic = False
                if len(nuke_alt.split(',')) == 2:
                    triallelic = True
                    altBP1, altBP2 = nuke_alt.split(',')
                elif len(nuke_alt.split(',')) == 3:
                    quadallelic = True
                    altBP1, altBP2, altBP3 = nuke_alt.split(',')

                if triallelic:
                    add.extend([chrome,pos,nuke_ref,altBP1, altBP2, np.nan]) # add values to list for output files
                elif quadallelic:
                    add.extend([chrome,pos,nuke_ref,altBP1, altBP2, altBP3])
                else:
                    add.extend([chrome,pos,nuke_ref, nuke_alt, np.nan, np.nan])

                for species in ids: #go species by species
                    alt = 0
                    ref = 0 #reset all variables after each species
                    tot = 0
                    missing = 0
                    alt2, alt3 = 0, 0

                    use = ids[species] #get list of ids
                    ## iterates through each species
                    for index in use: #go id by id for each species
                        data = lst[index][:3] #data for eact id
                        if len(data) < 2:
                            missing += 2
                        else:
                            ref += data.count('0')
                            alt += data.count('1')
                            missing += data.count('.')
                            if triallelic:
                                alt2 += data.count('2') #data for all species
                            if quadallelic:
                                alt2 += data.count('2')
                                alt3 +=  data.count('3')
                        tot += 2

                        if species in ['Pan','Pongo','Gorilla']: #puts in data depending on species (just ancestors)
                            if len(data) < 2:
                                ancest_missing += 2
                            else:
                                ancest_ref += data.count('0')
                                ancest_alt += data.count('1')
                                ancest_missing += data.count('.')
                                if triallelic:
                                    an_alt2 += data.count('2')
                                if quadallelic:
                                    an_alt2 += data.count('2')
                                    an_alt3 +=  data.count('3')
                            ancest_tot += 2


                    if len(use) != tot/2: #check 1
                        print("Ids don't match len of list")
                        exit()
                    if round(alt/tot+ref/tot+missing/tot+alt2/tot+alt3/tot,3) != 1: #check 2
                        print(chrome, pos, alt/tot+ref/tot+missing/tot+alt2/tot+alt3/tot)
                        print('Frequencies does not equal 1')
                        exit()

                    add.extend([ref/tot, alt/tot, alt2/tot, alt3/tot, missing/tot]) #add data to list



                ancest_missing_freq = ancest_missing/ancest_tot #get freq of missing in pan,pongo,gorilla
                freqs = [ancest_ref/ancest_tot, ancest_alt/ancest_tot, an_alt2/ancest_tot, an_alt3/ancest_tot, ancest_missing_freq]
                maxFreq = max(freqs)
                maxIndex = freqs.index(maxFreq)

                if ancest_missing_freq >= missing_stat or (maxIndex == 4):
                    ancest_state = 'NA' #when the ancestor nucleotide will be missing
                elif maxIndex == 0:
                    ancest_state = nuke_ref #when ancestor nucleotide will be the ref
                elif maxIndex == 1:
                    if triallelic or quadallelic:
                        ancest_state = altBP1
                    ancest_state = nuke_alt #when the ancestor nucleotide will be alt
                elif maxIndex == 2:
                    ancest_state = altBP2
                elif maxIndex == 3:
                    ancest_state = altBP3
                state_add.append(ancest_state)
                output_df_dic[len(output_df_dic)] = add #add list to dataframe
                state_df_dic[len(state_df_dic)] = state_add #add list to dataframe for ancestor state
    pd.DataFrame.from_dict(output_df_dic, orient='index', columns=cols).to_csv('/Frequencies/Frequencies.csv',index=False)  #convert dataframe to a text file
    pd.DataFrame.from_dict(state_df_dic,orient='index',columns=['Chromosome','Position','Nucleotide']).to_csv('/State/State.csv',index=False) #convert data frame to text file

makeTXT(vcf,Go,missing_stat)


