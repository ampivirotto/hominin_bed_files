import pandas as pd
import numpy as np
import sys
#import time

file =sys.argv[1]

def Ensembl(file):
    #start = time.time()
    #df = pd.DataFrame(columns=['Chromosome','Position','Ensembl'])
    dfdict = {}
    with open(file,'r') as data:
        for line in data:
            if line[0] == '>':
                lst = line.split(':')
                chrome = int(lst[2])
                pos= int(lst[3])
                end = int(lst[4])

            else:
                goat = line.strip('\n')
                for value in goat:
                    dfdict[pos] = [chrome, pos, value]
                    pos += 1
                    """
                    if pos%100000 == 0:
                        end = time.time()
                        print(pos, end-start)
                        start = end
                    """
        if pos-1 != end: #questionable on this test
            print(pos,end,sep='\t')
            print(int(end)-int(pos-1))

    pd.DataFrame.from_dict(dfdict, orient='index', columns=['Chromosome','Position','Ensembl']).to_csv('../ensembl/NoahEnsembl/'+chrome+'Ensembl.csv',index=False)
    #pd.DataFrame.from_dict(dfdict, orient='index', columns=['Chromosome','Position','Ensembl']).to_csv('./'+str(chrome)+'Ensembl.csv',index=False)


Ensembl(file)