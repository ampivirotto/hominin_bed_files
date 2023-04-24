import pandas as pd
import numpy as np
import sys 


file ='homo_sapiens_ancestor_1.fa' #sys.argv[1]

def Ensembl(file):
    df = pd.DataFrame(columns=['Chromosome','Position','Ensembl'])
    data = open(file,'r')
    for line in data:
        if line[0] == '>':
            lst = line.split(':')
            chrome = lst[2]
            pos= int(lst[3])
            end = int(lst[4])
            
        else:
            goat = str(line)
            for value in goat:
                df.loc[len(df.index)] = [chrome,pos,value]
                pos += 1
    if pos != end:
        print(pos,end,sep='\t')
        print(int(end)-int(pos))
        exit()
    df.to_csv(chrome+'Ensembl.txt',sep='\t',index=False)


Ensembl(file)