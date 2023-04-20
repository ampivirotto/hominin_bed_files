import pandas as pd
import numpy as np
import gzip

vcf = gzip.open('chr22_cds_primates.vcf.gz','r')
file = open('hominin_bed_files/primate_branch_ids.txt','r')
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
        species = list(df[df[0] == idName][1])[0]

        indDict[i] = species

        if species not in specDict:
            specDict[species] = [i]
        else:
            specDict[species].append(i)

    return indDict, specDict

for line in vcf:   #getting the ids add putting in a list
    line = line.decode('ASCII')
    if line[:3] == '#CH':
        no,Go = makeIndDictionary(line,file)
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



output_df = pd.DataFrame(columns=['Chromosome','Position','Ref','Alt','1000K-Ref','1000K-Alt','1000K-Missing','Pongo-Ref','Pongo-Alt','Pongo-Missing','Pan-Ref','Pan-ALt','Pan-Missing','Archaic-Ref','Archaic-Alt','Archaic-Missing'])
#above is a dataframe just to easily organize the data

totalLines = 77740
counter = 1 
for line in vcf:
    line = line.decode('ASCII')
    if round((counter / totalLines * 100),2) % 5 == 0:
        print(str(counter/totalLines * 100) + "%")
    counter += 1 
    if line[0] != '#':
        add = []
        lst = line.split()
        chrome = lst[0]
        pos = lst[1]
        nuke_ref = lst[3]
        nuke_alt = lst[4]
        add.extend([chrome,pos,nuke_ref,nuke_alt])
        for species in Go:
            alt = 0
            ref = 0 #reset all variables after each species
            tot = 0
            missing = 0
            use = Go[species] #get list of ids
            for index in use: #go id by id for each species
                data = lst[index][:3]
                if data == '0|0' or data == '0/0':
                    ref += 2
                    tot += 2
                elif data == '1|1' or data == '1/1': #record that data
                    alt += 2
                    tot += 2
                elif data == '0|1' or data == '1|0' or data == '0/1':
                    ref += 1
                    alt += 1
                    tot += 2
                else:
                    missing += 2
                    tot += 2
            if len(use) != tot/2:
                print("Ids don't match len of list")
                exit()
            if round(alt/tot+ref/tot+missing/tot,3) != 1:
                print(alt/tot+ref/tot+missing/tot)
                print('Frequencies does not equal 1')
                exit()
            add.extend([alt/tot,ref/tot,missing/tot]) #add data to list
    output_df.loc[len(output_df.index)] = add #add list to dataframe


output_df.to_csv('Freq.txt',sep='\t',index=False) #convert dataframe to a text file
