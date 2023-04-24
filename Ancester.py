import pandas as pd
import numpy as np
import gzip
import sys

missing_stat  = float(sys.argv[1]) #system input for missingness
vcf = sys.argv[2]
names = sys.argv[3]

vcf = gzip.open(vcf,'r')
idfile = open(names,'r')
# Id1000 = open('../ep/1000K.txt','r')
# IdPongo = open('../ep/Pongo.txt','r')
# IdPan = open('../ep/Pan.txt','r')
# IdArch = open('../ep/Archaics.txt','r')
"""
vcf = gzip.open('chr22_cds_primates.vcf.gz','r')
Id1000 = open('IDS/1000K.txt','r')
IdPongo = open('IDS/Pongo.txt','r')
IdPan = open('IDS/Pan.txt','r')
IdArch = open('IDS/Archaics.txt','r')
"""
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

# Go = {} #empty dictionary to get list of all different ids and species
# ktho = []
# for line in Id1000: #getting ids for 1000k into a dictionary
#     if len(line)>0:
#         Id = line.split()[0]
#         if Id in idlist:
#             ktho.append(idlist.index(Id))
# Go['1000K'] = ktho

# Arch = []
# for line in IdArch: #getting ids for archaics and putting into a dictionary
#     if len(line)>0:
#         Id = line.split()[0]
#         if Id in idlist:
#             Arch.append(idlist.index(Id))
# Go['Archaics'] = Arch

# Pan = []
# for line in IdPan: #getting ids fo pan and putting into dictionary
#     if len(line)>0:
#         Id = line.split()[0]
#         if Id in idlist:
#             Pan.append(idlist.index(Id))
# Go['Pan'] = Pan

# Pongo = []
# for line in IdPongo: #getting ids for pongo and putting into dictionary
#     line = line.strip()
#     if len(line)>0:
#         Id = line.split()[0]
#         if Id in idlist:
#             Pongo.append(idlist.index(Id))
# Go['Pongo'] = Pongo

def makeTXT(vcf,ids,missing):
    cols = ['Chromosome','Position','Ref','Alt', 'Alt2', 'Alt3']
    for species in ids.keys():
        for x in ['ref', 'alt', 'alt2', 'alt3', 'missing']:
            cols.append(species + "_" + x)
    output_df = pd.DataFrame(columns=cols)
    #above is a dataframe just to easily organize the data
    state_df = pd.DataFrame(columns=['Chromosome','Position','Nucleotide']) #Ensembl will be added later

    # totalLines = 77740
    # counter = 1
    for line in vcf:
        line = line.decode('ASCII')
        # if round((counter / totalLines * 100),2) % 5 == 0:
        #     print(str(counter/totalLines * 100) + "%")
        # counter += 1
        if line[0] != '#':
            add = [] #lest to add to freq file
            state_add = [] #list to add to state file
            lst = line.split() #list of values in line
            chrome = lst[0] #chrom in line
            pos = lst[1] #position in line
            nuke_ref = lst[3] #ref in line
            nuke_alt = lst[4] #alt in line
            ancest_missing, ancest_tot, ancest_alt, ancest_ref = 0, 0, 0, 0
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
                            alt2 += data.count('2')
                        if quadallelic:
                            alt2 += data.count('2')
                            alt3 +=  data.count('3')
                    tot += 2

                    if species in ['Pan','Pongo','Gorilla']: #puts in data depeending on species
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

                    """
                    if data == '0|0' or data == '0/0':
                        ref += 2
                        tot += 2
                        if species in ['Pan','Pongo','Gorilla']: #puts in data depeending on species
                            ancest_ref += 2
                            ancest_tot += 2
                    elif data == '1|1' or data == '1/1': #record that data
                        alt += 2
                        tot += 2
                        if species in ['Pan','Pongo','Gorilla']:
                            ancest_alt += 2
                            ancest_tot += 2
                    elif data == '0|1' or data == '1|0' or data == '0/1':
                        ref += 1
                        alt += 1
                        tot += 2
                    else:
                        missing += 2
                        tot += 2
                        if species in ['Pan','Pongo','Gorilla']:
                            ancest_missing += 2
                            ancest_tot += 2
                    """
                if len(use) != tot/2: #check 1
                    print("Ids don't match len of list")
                    exit()
                if round(alt/tot+ref/tot+missing/tot,3) != 1: #check 2
                    print(alt/tot+ref/tot+missing/tot)
                    print('Frequencies does not equal 1')
                    exit()
                if not triallelic:
                    alt2 = 0
                if not quadallelic:
                    alt3 = 0
                add.extend([ref/tot, alt/tot, alt2/tot, alt3/tot, missing/tot]) #add data to list

            if not triallelic:
                an_alt2 = 0
            if not quadallelic:
                an_alt3 = 0
            ancest_missing_freq = ancest_missing/ancest_tot #get freq of missing in pan,pongo,gorilla
            #ancest_alt_freq = ancest_alt/ancest_tot #get freq of alt in pan,pongo,gorilla
            #ancest_ref_freq = ancest_ref/ancest_tot #freq of ref in pan,pongo,gorilla
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
            output_df.loc[len(output_df.index)] = add #add list to dataframe
            state_df.loc[len(state_df.index)] = state_add #add list to dataframe for ancestor state
    output_df.to_csv('Freq.txt',sep='\t',index=False) #convert dataframe to a text file
    state_df.to_csv('State.txt',sep='\t',index=False) #convert data frame to text file

makeTXT(vcf,Go,missing_stat)


# output_df = pd.DataFrame(columns=['Chromosome','Position','Ref','Alt','1000K-Ref','1000K-Alt','1000K-Missing','Pongo-Ref','Pongo-Alt','Pongo-Missing','Pan-Ref','Pan-ALt','Pan-Missing','Archaic-Ref','Archaic-Alt','Archaic-Missing'])
# #above is a dataframe just to easily organize the data
# state_df = pd.DataFrame(columns=['Chromosome','Position','Nucleotide']) #Ensembl will be added later

# totalLines = 77740
# counter = 1
# for line in vcf:
#     line = line.decode('ASCII')
#     if round((counter / totalLines * 100),2) % 5 == 0:
#         print(str(counter/totalLines * 100) + "%")
#     counter += 1
#     if line[0] != '#':
#         add = []
#         state_add = []
#         lst = line.split()
#         chrome = lst[0]
#         pos = lst[1]
#         nuke_ref = lst[3]
#         nuke_alt = lst[4]
#         ancest_missing = 0
#         ancest_tot = 0
#         ancest_alt = 0
#         ancest_ref = 0
#         add.extend([chrome,pos,nuke_ref,nuke_alt])
#         state_add.extend([chrome,pos])
#         for species in Go:
#             alt = 0
#             ref = 0 #reset all variables after each species
#             tot = 0
#             missing = 0
#             use = Go[species] #get list of ids
#             for index in use: #go id by id for each species
#                 data = lst[index][:3]
#                 if data == '0|0' or data == '0/0':
#                     ref += 2
#                     tot += 2
#                     if species in ['Pan','Pongo','Gorilla']:
#                         ancest_ref += 2
#                         ancest_tot += 2
#                 elif data == '1|1' or data == '1/1': #record that data
#                     alt += 2
#                     tot += 2
#                     if species in ['Pan','Pongo','Gorilla']:
#                         ancest_alt += 2
#                         ancest_tot += 2
#                 elif data == '0|1' or data == '1|0' or data == '0/1':
#                     ref += 1
#                     alt += 1
#                     tot += 2
#                 else:
#                     missing += 2
#                     tot += 2
#                     if species in ['Pan','Pongo','Gorilla']:
#                         ancest_missing += 2
#                         ancest_tot += 2
#             if len(use) != tot/2: #check 1
#                 print("Ids don't match len of list")
#                 exit()
#             if round(alt/tot+ref/tot+missing/tot,3) != 1: #check 2
#                 print(alt/tot+ref/tot+missing/tot)
#                 print('Frequencies does not equal 1')
#                 exit()
#             add.extend([alt/tot,ref/tot,missing/tot]) #add data to list
#         ancest_missing_freq = ancest_missing/ancest_tot
#         ancest_alt_freq = ancest_alt/ancest_tot
#         ancest_ref_freq = ancest_ref/ancest_tot
#         if ancest_missing_freq >= missing_stat or (ancest_missing_freq>ancest_alt_freq and ancest_missing_freq > ancest_ref_freq):
#             ancest_state = 'NA'
#         elif ancest_ref_freq > ancest_missing_freq and ancest_ref_freq > ancest_alt_freq:
#             ancest_state = nuke_ref
#         elif ancest_alt_freq > ancest_missing_freq and ancest_alt_freq > ancest_ref_freq:
#             ancest_state = nuke_ref
#         state_add.append(ancest_state)
#         output_df.loc[len(output_df.index)] = add #add list to dataframe
#         state_df.loc[len(state_df.index)] = state_add

# #output of textfile for chromosome, position, nucleotide, ensembl (European ancester state)
# #calc most common allele for pongo, pan, gorilla /also get the missing amount

# output_df.to_csv('Freq.txt',sep='\t',index=False) #convert dataframe to a text file

# state_df.to_csv('State.txt',sep='\t',index=False)
