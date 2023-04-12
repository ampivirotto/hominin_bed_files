from math import *
import numpy as np
import pandas as pd



file = 'refGene.txt' #'/Users/noahpeles/Downloads/refGene.txt'
file2 = 'refseq_ids_w_eps.csv' #'/Users/noahpeles/Coding/Reseach/hominin_bed_files/refseq_ids_w_eps.csv'

## goes through file and takes all refseq ids and adds to a list
ToUse = []
for line in open(file2):
    ls = line.strip()
    ToUse.append(ls)

identity = []
chrom = []
starts = []
ends = []
strands = []

## goes through file of genes and makes bed file line if the id is in the list
for line in open(file,'r'):
    ls = line.strip().split()
    if ls[1] in ToUse:
        begins = ls[9].split(',')
        stops = ls[10].split(',')
        for i in range(len(begins)-1):
            identity.append(ls[1])
            chrom.append(ls[2])
            ends.append(stops[i])
            strands.append(ls[3])
            if ls[3]=='-':
                starts.append(int(begins[i])-1)
            else:
                starts.append(begins[i])
    else:
        continue

chom = []
for k in chrom:
    ToCheck = k[3:5].replace('_','')
    if ToCheck.isalpha() == False:
        chom.append(int(ToCheck))
    else:
        chom.append(ToCheck)


df_Both = pd.DataFrame({'Chromosome':chom,'Identity':identity,'Strand':strands,'Start Points':starts,'End Points':ends})
#all = len(df_Both.index) #To check data
print(df_Both)
print(len(df_Both.index))
df_dropped = df_Both.drop_duplicates(subset=['Chromosome','Start Points','End Points']) #table with dropped duplicates
df_dropped = df_dropped.sort_values(by=['Chromosome','Start Points'])
df_dropped = df_dropped.drop(columns=['Identity','Strand'])

#chrm22 = df_dropped.loc[df_dropped['Chromosome']==22]
#chrm22.to_csv('22.bed',sep='\t',header=False,index=False)

#new_chrom = []
#for k in df_dropped['Chromosome']:
    #good = str(k)
    #new_chrom.append('chr'+good)

#df_dropped['Chromosome'] = new_chrom
#df_dropped.to_csv('EP_CDS.bed',sep='\t',index=False,header=False)
#df_dropped.to_csv('EP_CDS.csv',index=False)

#df_d = df_Both[df_Both.duplicated(subset=['Chromosome','Start Points','End Points'],keep=False)] #Table of what the duplicates were
#df_d = df_d.sort_values(by=['Chromosome','Start Points'])
#df_d = df_d.drop(columns=['Strand'])
#df_d.to_csv('Duplicated.bed',sep='\t',index=False)

#dropped1 = len(df_d['Start Points'].unique())
#dropped2 = len(df_d.index)
#print(all)                     #To check data numbers
#print(len(df_dropped.index))
#print(dropped2-dropped1)

#df_Both = df_Both.drop(columns=['Identity','Strand']) #Table of everything
#df_Both = df_Both.sort_values(by=['Chromosome','Start Points'])
#df_Both.to_csv('EP_CDS_All.bed',sep='\t',index=False)


for i in df_dropped['Chromosome'].unique():    #Making table for each chomosome
    df = df_dropped.loc[df_dropped['Chromosome']==i]
    print(len(df.index))
    name = str(i)+'.bed'
    df.to_csv("../ChromosomeFiles/" + name,sep='\t',index=False,header=False)