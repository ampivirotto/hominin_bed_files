import pandas as pd
import random


New = open('open.vcf',"w")
NonCall = pd.DataFrame({'Chr':[22,22,22],'Start':[49,102,176],'End':[77,104,200]})
dic = {}
for i in range(len(NonCall.index)):
    row = NonCall.iloc[i]
    for k in range(row[1],row[2]+1):
        dic[k] = random.choice(['A','T','G','C'])

addin_df = pd.DataFrame(columns=['Chrom','Pos','rsid','Ref','Alt','Qual','Info','Indi'])
for i in range(len(NonCall.index)):
    row = NonCall.iloc[i]
    for k in range(row[1],row[2]+1):
        ref = dic[k]
        go = [row[0],k,'-',ref,'N','-','-','./.']
        addin_df.loc[len(addin_df.index)] = go


start = NonCall.iloc[0][1]
get = 0
with open('/Users/noahpeles/Downloads/VCFTool.vcf', 'r') as file:
    for line in file:
        if line[0] != '#':
            LineParts = line.strip().split('\t')
            if int(LineParts[1]) > start:
                for i in range(start,NonCall.iloc[get][2]):
                    new_line = '    '.join(addin_df.loc[addin_df['Pos']==i])
                    New.write(new_line)
                start = NonCall.iloc[get][1]
                get += 1
            else:
                New.write(line)

