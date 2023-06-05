import sys
import pysam

chromosome = sys.argv[1] #only input will be the number of a chomosome


State_path = './States/chrome_{}.csv'.format(chromosome)
bedfile = './bedfile/{}.merged.bed'.format(chromosome)


def StateDictionary(file2Path):
    State = {}
    with open(file2Path,'r') as data:
        for line in data:
            lis = line.strip('\n').split(',')
            if lis[1][0].isnumeric():
                Position = int(lis[1])
                State[Position] = lis[2]
    return State

def REF(chrnum,bedfile):
    """
    get reference genome
    """
    refDict = {} 
    genome = pysam.FastaFile('./hg19_reference/chr{}.fa.gz'.format(chrnum)) #reading in file using pysam
    with open(bedfile) as f:
        for line in f:
            splitLine = line.strip('\n').split('\t')
            for x in range(int(splitLine[1]), int(splitLine[2])+1):
                pos = genome.fetch(genome.references[0], start=x, end=x+1) #fetching the position
                refDict[x] = pos
    return refDict



def MakeVCF(bed,Spath,chromosome):
    '''Make VCF file for the state and ensembl data, take in the bed file for iteration
    and paths for the various different files'''

    Sdic = StateDictionary(Spath)
    Rdic = REF(chromosome,bed)

    output = open('Ancestor_chr{}.vcf'.format(chromosome),'w')
    output.write('##fileformat=VCFv4.2'+'\n')
    output.write('##Database=PrimateIndividuals'+'\n')
    headers2 = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','PRIMATE'])
    output.write(headers2+'\n')


    with open(bed,'r') as file: #could interate through the merged bed file, and go through those ranges
        for line in file:
                lis = line.strip('\n').split('\t')
                start,end = int(lis[1]),int(lis[2])
                for position in range(start,end+1):
                    Ancestor_Add = [chromosome]
                    Ref = Rdic[position]
                    Ancestor_Add += [str(position),'.',Ref]
                    Ancestor_nuke = Sdic.get(position,'NA')
                    if Ancestor_nuke == Ref:
                        Ancestor_Add += ['.','.','.','.','.','0|0']
                    elif Ancestor_nuke == 'NA':
                       Ancestor_Add += ['.','.','.','.','.','.|.'] 
                    else:
                        if Ancestor_nuke.isalpha():
                            Ancestor_Add += [Ancestor_nuke,'.','.','.','.','1|1']
                        elif Ancestor_nuke == '.':
                            Ancestor_Add += ['.','.','.','.','.','.']
                    Ancestor_String = '\t'.join(Ancestor_Add)
                    output.write(Ancestor_String+'\n')
    output.close()

MakeVCF(bedfile,State_path,chromosome)