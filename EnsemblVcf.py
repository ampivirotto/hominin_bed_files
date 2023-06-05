import sys
import pysam


chromosome = sys.argv[1] #only input will be the number of a chomosome

Ensembl_path = './../ensembl/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{}.fa'.format(chromosome)
bedfile = './bedfiles/{}.merged.bed'.format(chromosome)

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
                pos = genome.fetch(genome.references[0], start=x-1, end=x) #fetching the position   ## AP fix bug in finding reference (off by one error)
                refDict[x] = pos.upper()  ## AP make sure it's uppercase
    return refDict



def MakeVCF(bed,Epath,chromosome):
    '''Make VCF file for the state and ensembl data, take in the bed file for iteration
    and paths for the various different files'''

    Edic = Ensembl(Epath)
    Rdic = REF(chromosome,bed)
    output1 = open('EnsemblVcf/Ensembl{}.vcf'.format(chromosome),'w')
    output1.write('##fileformat=VCFv4.2'+'\n')
    output1.write('##Database=Ensembl'+'\n')
    headers1 = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','ENSEMBL'])
    output1.write(headers1+'\n')


    with open(bed,'r') as file: #could interate through the merged bed file, and go through those ranges
        for line in file:
                lis = line.strip('\n').split('\t')
                start,end = int(lis[1]),int(lis[2])
                for position in range(start,end+1):
                    Ensembl_Add = [chromosome]
                    Ref = Rdic[position]
                    Ensembl_Add += [str(position),'.',Ref]
                    Ensembl_nuke = Edic.get(position,'NA').upper()
                    if Ensembl_nuke == Ref:
                        Ensembl_Add += ['.','.','.','.','.','0|0']
                    elif Ensembl_nuke == 'NA':
                        Ensembl_Add += ['.','.','.','.','.','.|.']
                    else:
                        if Ensembl_nuke.isalpha():
                            Ensembl_Add += [Ensembl_nuke,'.','.','.','.','1|1']  ## AP make it uppercase
                        elif Ensembl_nuke == '.':
                            Ensembl_Add += ['.','.','.','.','.','.']
                    Ensembl_String = '\t'.join(Ensembl_Add)
                    output1.write(Ensembl_String+'\n')
    output1.close()

MakeVCF(bedfile,Ensembl_path,chromosome)