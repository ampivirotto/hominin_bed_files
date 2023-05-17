import sys

chromosome = sys.argv[1] #only input will be the number of a chomosome

Ensembl_path = '../ensembl/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{}.fa'.format(chromosome)
State_path = '/States/chrome_{}.csv'.format(chromosome) #made up name here
original_vcf = 'Merged_chr{}.csv'.format(chromosome) #made up name here

def Ensembl(filePath):
    dfdict = {}
    with open(filePath,'r') as data:
        for line in data:
            if line[0] == '>':
                lst = line.split(':')
                pos= int(lst[3])
                end = int(lst[4])

            else:
                goat = line.strip('\n')
                for value in goat:
                    dfdict[pos] = [value]
                    pos += 1
        if pos-1 != end: #questionable on this test
            print(pos,end,sep='\t')
            print(int(end)-int(pos-1))
    return dfdict

def StateDictionary(file2Path):
    State = {}
    with open(file2Path,'r') as data:
        for line in data:
            lis = line.split(',')
            if lis[1][0].isnumeric():
                Position = lis[1]
                State[Position] = lis[2]
    return State

EnsemblDic = Ensembl(Ensembl_path)
StateDic = StateDictionary(State_path)

output1 = open('/StateVcf/Ensembl{}.vcf'.format(chromosome),'w')
output1.write('##fileformat=VCFv4.2'+'\n')
output1.write('##Database=Ensembl'+'\n')
headers1 = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','ENSEMBL'])
output1.write(headers1+'\n')

output2 = open('/StateVcf/Ancestor_chr{}.vcf'.format(chromosome),'w')
output2.write('##fileformat=VCFv4.2'+'\n')
output2.write('##Database=PrimateIndividuals'+'\n')
headers2 = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','PRIMATE'])
output2.write(headers2+'\n')

with open(original_vcf,'r') as file: #could interate through the merged bed file, and go through those ranges
    for line in file: 
        if line[0] != '#':
            lis = line.strip().split()
            Ensembl_Add,Ancestor_Add = [chromosome],[chromosome]
            position,Ref = lis[1],lis[3]
            Ensembl_Add += [position,'.',Ref]
            Ancestor_Add += [position,'.',Ref]
            Ensembl_nuke,Ancestor_nuke = EnsemblDic[position],StateDic[position]
            if Ensembl_nuke == Ref:
                Ensembl_Add += ['.','.','.','.','.','0|0']
            else:
                if Ensembl_nuke.isalpha():
                    Ensembl_Add += [Ensembl_nuke,'.','.','.','.','1|1']
                elif Ensembl_nuke == '.':
                    Ensembl_Add += ['.','.','.','.','.','.']
            if Ancestor_nuke == Ref:
                Ancestor_Add += ['.','.','.','.','.','0|0']
            else:
                if Ancestor_nuke.isalpha():
                    Ancestor_Add += [Ensembl_nuke,'.','.','.','.','1|1']
                elif Ancestor_nuke == '.':
                    Ancestor_Add += ['.','.','.','.','.','.']
            Ensembl_String = '\t'.join(Ensembl_Add)
            Ancestor_String = '\t'.join(Ancestor_Add)
            output1.write(Ensembl_String+'\n')
            output2.write(Ancestor_String+'\n')

#look into switching to merged file instead iterating through vcf