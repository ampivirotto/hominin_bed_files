import pandas as pd
import numpy as np
import pickle



New = open('chr22.vcf',"w")

FirstDF = pd.read_csv('/Users/noahpeles/Coding/Reseach/AddinFiles/chr22_Pongo_abelli.bed',sep='\t', header = None)
#print(FirstDF)
#FirstDF.columns = ['CHROM','Start','End','?']

NonCall = FirstDF#.loc[FirstDF['CHROM']=='chr22']

dic = pd.read_pickle(r'/Users/noahpeles/Coding/Reseach/AddinFiles/ensemble_grch37_ancestor_22.p')

missing = open('missingSites.txt','w')
addin_df = pd.DataFrame.from_dict(dic, orient='index').reset_index()
print(addin_df)
"""
addin_df = pd.DataFrame(columns=['CHROM','POS','ID','REF','ALT','QUAL', 'FILTER', 'INFO', 'FORMAT','Indi'])
for i in range(len(NonCall.index)):
    row = NonCall.iloc[i]
    for k in range(row[1],row[2]+1):
        try: 
            ref = dic[k]
            go = [str(row[0]),str(k),'.',ref,'.','.','.', '.', '.','./.']
            addin_df.loc[len(addin_df.index)] = go
        except:
            print(k)
            go = [str(row[0]),str(k),'.','.','.','.','.','.','.','./.']
            addin_df.loc[len(addin_df.index)] = go
            missing.write(str(k)+'\n')
missing.close()
"""

start = int(NonCall.iloc[0][1])
go_to = int(NonCall.iloc[0][2])+1
get = 0 
end = NonCall.iloc[len(NonCall)-1][2]
with open('/Users/noahpeles/Coding/Reseach/AddinFiles/chr22_cds_Pongo_abelii_filtered.vcf', 'r') as file:
    for line in file:
        if line[0] != '#':
            LineParts = line.strip().split('\t')
            if int(LineParts[1]) > start:
                for i in range(start,go_to):
                    new_line = '\t'.join(addin_df.loc[addin_df['POS'].astype(int) == i].values.flatten().tolist())
                    New.write(new_line +'\n')
                #New.write(line)
                if int(LineParts[1]) < end: #avoids errors for the last add in
                    get += 1
                    start = int(NonCall.iloc[get][1])
                    go_to = int(NonCall.iloc[get][2])+1
            else:
                New.write(line)

New.close()