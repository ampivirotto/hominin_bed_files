

import pandas as pd
import numpy as np
import gzip
import pickle
import sys


def turnBedintoDF(bed):
    bedDict = {}
    count = 0
    for i in range(len(bed.index)):
        row = bed.iloc[i]
        for pos in np.arange(row[1], row[2]):
            bedDict[count] = [row[0], pos]
            count += 1
    return pd.DataFrame.from_dict(bedDict, orient = 'index')

def addVCFcols(df, newCol):
    df.loc[:, newCol] = ['.' for i in range(len(df))]
    return df

def addVCFindividuals(df, numInds):
    for x in range(numInds):
        df.loc[:, 'ind' + str(x)] = ['./.' for i in range(len(df))]
    return df

def main(chromNum, ncBED, vcffname):
    ## create new vcf output file
    outfname = vcffname.strip('.vcf.gz') + '_nonCallRegionAdded.vcf'
    New = open(outfname,"w")

    ## read in bed file of noncallable regions
    FirstDF = pd.read_csv(ncBED, sep='\t', header = None)
    NonCall = turnBedintoDF(FirstDF)

    ## read in dictionary of ancestor
    dic = pd.read_pickle(r'grch37_poskey_baseval_dic_{}_based.p'.format(chromNum))

    ## make file for writing out missing sites
    #missing = open('missingSites.txt','w')

    ## turn dictionary into dataframe
    ensemble_df = pd.DataFrame.from_dict(dic, orient='index').reset_index()

    ## select rows for noncallable
    addin_df = ensemble_df.merge(NonCall, left_on = 'index', right_on = 1)

    ## turn dfs into a row df for each site in vcf format
    addin_df = addin_df.rename(columns = {'0_y':'CHROM', 'index':'POS', '0_x':'REF'})
    for col in ['ID','ALT','QUAL', 'FILTER', 'INFO', 'FORMAT']:
        addin_df = addVCFcols(addin_df, col)
    addin_df['CHROM'] = addin_df['CHROM'].str.strip('chr')
    addin_df = addin_df[['CHROM','POS','ID','REF','ALT','QUAL', 'FILTER', 'INFO', 'FORMAT']]
    addin_df = addVCFindividuals(addin_df, numInds)
    addin_df = addin_df.astype('str', copy = 'False')
    print(len(addin_df))

    addin_df = addin_df.drop_duplicates().reset_index(drop = True)
    #running_index = 0
    #lastline = max(FirstDF[1])+ 1
    with gzip.open(vcffname, 'r') as file:
        for line in file:
            line = line.decode('ASCII')
            if line[0] != '#':
                addin_df.loc[len(addin_df.index)] = line.strip('\n').split("\t")
            else:
                New.write(line)

    addin_df.drop_duplicates(subset = ['CHROM', 'POS'], keep = 'last', inplace = True)
    addin_df.sort_values(by = ['CHROM', 'POS'], inplace = True)
    addin_df = addin_df.reset_index(drop = True)

    print(len(addin_df))


    # --- Alexander was here ---
    for row in range(len(addin_df)):
        new_line = '\t'.join(addin_df.loc[row].values.flatten().tolist())
        New.write(new_line + '\n')

    New.close()

if __name__ == '__main__':
    ## python NAME.py CHROM BEDFILE VCFFILE
    main(sys.argv[1], sys.argv[2], sys.argv[3])

## CODE NOT USED
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

"""
LineParts = line.strip().split('\t')
            pos = int(LineParts[1])
            dflist = []
            for x in np.arange(running_index, len(FirstDF)):
                start = int(FirstDF.iloc[x][1])
                go_to = int(FirstDF.iloc[x][2])+1
                if (start <= 24179249) & (go_to > 24179249):
                    print('here')
                #if (pos == start) or (pos == go_to):
                    #print(pos)
                if pos > go_to:
                    if lastline > start:
                        dflist.append(addin_df[(addin_df['POS'].astype(int) > lastline) & (addin_df['POS'].astype(int) <= go_to)])
                        lastline = go_to
                    else:
                        dflist.append(addin_df[(addin_df['POS'].astype(int) >= start) & (addin_df['POS'].astype(int) <= go_to)])
                        lastline = go_to
                elif (pos > start) & (pos < go_to):
                    if pos > lastline:
                        dflist.append(addin_df[(addin_df['POS'].astype(int) > lastline) & (addin_df['POS'].astype(int) < pos)])
                    else:
                        dflist.append(addin_df[(addin_df['POS'].astype(int) >= start) & (addin_df['POS'].astype(int) < pos)])
                    nonCallDF = pd.concat(dflist).drop_duplicates(subset='POS').sort_values(by = 'POS').reset_index(drop = True)
                    for j in range(len(nonCallDF)):
                        new_line = '\t'.join(nonCallDF.loc[j].values.flatten().tolist())
                        if nonCallDF.loc[j]['POS'] == 17443044:
                            print('here')
                        New.write(new_line +'\n')
                    New.write(line)
                    lastline = pos
                    running_index = x
                    break
                elif (pos < start):
                    if len(dflist) > 0:
                        nonCallDF = pd.concat(dflist).reset_index(drop = True).sort_values(by = 'POS')
                        for j in range(len(nonCallDF)):
                            new_line = '\t'.join(nonCallDF.loc[j].values.flatten().tolist())
                            if nonCallDF.loc[j]['POS'] == 17443044:
                                print('here')
                            New.write(new_line +'\n')
                    New.write(line)
                    lastline = pos
                    running_index = x
                    break
"""