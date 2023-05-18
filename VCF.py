import sys
import pysam

chromosome = sys.argv[1] #only input will be the number of a chomosome

Ensembl_path = './../ensembl/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{}.fa'.format(chromosome)
State_path = './States/chrome_{}.csv'.format(chromosome)   
bedfile = './bedfile/{}.merged.bed'.format(chromosome)

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
                    dfdict[pos] = value
                    pos += 1
        if pos-1 != end:
            print(pos,end,sep='\t')
            print(int(end)-int(pos-1))
    return dfdict

def StateDictionary(file2Path):
    State = {}
    with open(file2Path,'r') as data:
        for line in data:
            lis = line.strip('\n').split(',')
            if lis[1][0].isnumeric():
                Position = int(lis[1])
                State[Position] = lis[2]
    return State

def REF(chrnum):
    """
    get reference genome
    """
    refDict = {}
    genome = pysam.FastaFile('./hg19_reference/chr{}.fa.gz ./bedfiles/{}.merged.bed'.format(chrnum)) #reading in file using pysam
    with open('hominin_bed_files/output/{}.merged.bed'.format(chrnum)) as f:
        for line in f:
            splitLine = line.strip('\n').split('\t')
            for x in range(int(splitLine[1]), int(splitLine[2])+1):
                pos = genome.fetch(genome.references[0], start=x, end=x+1) #fetching the position
                refDict[x] = pos
    return refDict



def MakeVCF(bed,Epath,Spath,chromosome):
    '''Make VCF file for the state and ensembl data, take in the bed file for iteration
    and paths for the various different files'''

    Edic = Ensembl(Epath)
    Sdic = StateDictionary(Spath)
    Rdic = REF(chromosome)
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


    with open(bed,'r') as file: #could interate through the merged bed file, and go through those ranges
        for line in file:
                lis = line.strip('\n').split('\t')
                start,end = int(lis[1]),int(lis[2])
                for position in range(start,end+1):
                    Ensembl_Add,Ancestor_Add = [chromosome],[chromosome]
                    Ref = Rdic[position]
                    Ensembl_Add += [str(position),'.',Ref]
                    Ancestor_Add += [str(position),'.',Ref]
                    Ensembl_nuke,Ancestor_nuke = Edic.get(position,'NA'),Sdic.get(position,'NA')
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
                            Ancestor_Add += [Ancestor_nuke,'.','.','.','.','1|1']
                        elif Ancestor_nuke == '.':
                            Ancestor_Add += ['.','.','.','.','.','.']
                    Ensembl_String = '\t'.join(Ensembl_Add)
                    Ancestor_String = '\t'.join(Ancestor_Add)
                    output1.write(Ensembl_String+'\n')
                    output2.write(Ancestor_String+'\n')
    output1.close()
    output2.close()

MakeVCF(bedfile,Ensembl_path,State_path,chromosome)
