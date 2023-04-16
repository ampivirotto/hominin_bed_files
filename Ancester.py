import pandas as pd
import numpy as np
import gzip

vcf = gzip.open('../chr22/chr22_cds_primates.vcf.gz','r')
Id1000 = open('../ep/1000K.txt','r')
IdPongo = open('../ep/Pongo.txt','r')
IdPan = open('../ep/Pan.txt','r')
IdArch = open('../ep/Archaics.txt','r')
"""
vcf = gzip.open('chr22_cds_primates.vcf.gz','r')
Id1000 = open('IDS/1000K.txt','r')
IdPongo = open('IDS/Pongo.txt','r')
IdPan = open('IDS/Pan.txt','r')
IdArch = open('IDS/Archaics.txt','r')
"""

for line in vcf:
    line = line.decode('ASCII')
    if line[:3] == '#CH':
        idlist = line.split()
        break

Go = {}
ktho = []
for line in Id1000:
    if len(line)>0:
        Id = line.split()[0]
        if Id in idlist:
            ktho.append(idlist.index(Id))
Go['1000K'] = ktho

Arch = []
for line in IdArch:
    if len(line)>0:
        Id = line.split()[0]
        if Id in idlist:
            Arch.append(idlist.index(Id))
Go['Archaics'] = Arch

Pan = []
for line in IdPan:
    if len(line)>0:
        Id = line.split()[0]
        if Id in idlist:
            Pan.append(idlist.index(Id))
Go['Pan'] = Pan

Pongo = []
for line in IdPongo:
    line = line.strip()
    if len(line)>0:
        Id = line.split()[0]
        if Id in idlist:
            Pongo.append(idlist.index(Id))
Go['Pongo'] = Pongo



output_df = pd.DataFrame(columns=['Chromosome','Position','Ref','Alt','1000K-Ref','1000K-Alt','1000K-Missing','Pongo-Ref','Pongo-Alt','Pongo-Missing','Pan-Ref','Pan-ALt','Pan-Missing','Archaic-Ref','Archaic-Alt','Archaic-Missing'])

totalLines = 546653
counter = 1
for line in vcf:
    line = line.decode('ASCII')
    if (counter / totalLines * 100) % 5 == 0.0:
        print(str(round(counter/totalLines * 100, 3)) + "%")
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
            ref = 0
            tot = 0
            missing = 0
            use = Go[species]
            for index in use:
                data = lst[index][:3]
                if data == '0|0' or data == '0/0':
                    ref += 2
                    tot += 2
                elif data == '1|1' or data == '1/1':
                    alt += 2
                    tot += 2
                elif data == '0|1' or data == '1|0' or data == '0/1':
                    ref += 1
                    alt += 1
                    tot += 2
                else:
                    missing += 2
                    tot += 2
            add.extend([alt/tot,ref/tot,missing/tot])
    output_df.loc[len(output_df.index)] = add


output_df.to_csv('Freq.txt',sep='\t',index=False)
