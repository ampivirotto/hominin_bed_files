import pandas as pd
import numpy as np
import random


New = open('open.vcf',"w")
NonCall = pd.DataFrame({'Chr':[22,22,22],'Start':[49,102,176],'End':[77,104,200]})
dic = {}
for i in range(len(NonCall.index)):
    row = NonCall.iloc[i]
    for k in range(row[1],row[2]+1):
        dic[k] = random.choice(['A','T','G','C'])

addin_df = pd.DataFrame(columns=['CHROM','POS','ID','REF','ALT','QUAL', 'FILTER', 'INFO', 'FORMAT','Indi'])
for i in range(len(NonCall.index)):
    row = NonCall.iloc[i]
    for k in range(row[1],row[2]+1):
        ref = dic[k]
        go = [str(row[0]),str(k),'.',ref,'.','.','.', '.', '.','./.']
        addin_df.loc[len(addin_df.index)] = go


start = NonCall.iloc[0][1]
get = 0
end = NonCall.iloc[len(NonCall)-1][2]
with open('VCFTool.vcf', 'r') as file:
    for line in file:
        if line[0] != '#':
            LineParts = line.strip().split('\t')
            if int(LineParts[1]) > start and int(LineParts[1]) < end:
                try:
                    for i in range(start,NonCall.iloc[get][2]):
                        new_line = '\t'.join(addin_df.loc[addin_df['POS'].astype(int) == i].values.flatten().tolist())
                        New.write(new_line + '\n')
                    start = NonCall.iloc[get][1]
                    get += 1
                except:
                    print('error')
            else:
                New.write(line)

New.close()
