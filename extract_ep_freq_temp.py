#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     11/04/2023
# Copyright:   (c) tuk32868 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import gzip
import pandas as pd
import numpy as np
import pysam
import Bio.Seq as seq
import sys


def readVCF(vcffn, indLoc, o, ep, iso, fa):
    with open(vcffn) as f:
        for line in f:
            if not line.startswith('#'):
                readLine(line, indDict, spD, o, ep, iso, fa)
            elif line.startswith('#C'):
                indDict, spD = makeIndDictionary(line, indLoc)


def readGZIPvcf(vcffn, indLoc, o, ep, iso, fa):
    with gzip.open(vcffn) as f:
        for line in f:
            line = line.decode('ASCII')
            if not line.startswith("#"):
                readLine(line, indDict, spD, o, ep, iso, fa)
            elif line.startswith('#CHROM'):
                indDict, spD = makeIndDictionary(line, indLoc)


def makeIndDictionary(line, location):
    indDict = {}
    specDict = {}

    df = pd.read_csv(location + 'primate_branch_ids.txt', sep = '\t', header = None)

    lineSplit = line.rstrip('\n').split('\t')

    for i in range(9, len(lineSplit)):
        idName = lineSplit[i]
        species = list(df[df[0] == idName][1])[0]

        indDict[i] = species

        if species not in specDict:
            specDict[species] = [i]
        else:
            specDict[species].append(i)

    return indDict, specDict


def checkFreq(ids, line, secondAlt = False):
    missing, ref, alt = 0, 0, 0
    if secondAlt:
        alt2 = 0

    for idi in ids:
        gt = line[idi][:3]

        if len(gt) == 1:
            missing += 2
        else:
            alt += gt.count('1')
            ref += gt.count('0')
            missing = gt.count('.')
            if secondAlt:
                alt2 = gt.count('2')
    totalGenomes = len(ids) * 2
    if secondAlt:
        return [ref/totalGenomes, alt/totalGenomes, missing/totalGenomes, alt2/totalGenomes]
    return [ref/totalGenomes, alt/totalGenomes, missing/totalGenomes]

def checkstrand(tbl):
    positions = list(tbl['nuc_pos'])
    if positions[0] > positions[1]:
        assert positions[0] - 1 == positions[1]
        return '-'
    elif positions[0] < positions[1]:
        assert positions[0] + 1 == positions[1]
        return '+'
    else:
        print('ERROR STRAND NOT IDENTIFIED')
        exit()

def readLine(line, indDict, specs, o, ep, iso, fa):
    splitLine = line.strip('\n').split('\t')

    lineDict = {}

    pos = splitLine[1]
    chrnum = splitLine[0]
    ref = splitLine[3]
    alt = splitLine[4]

    if 'chr' in chrnum:
        isochr = chrnum
    else:
        isochr = 'chr' + str(chrnum)

    if len(alt.split(',')) > 2:
        return
    elif (len(alt.split(',')) == 2) or (len(alt.split(','))== 1):
        outProps = checkFreq(specs['mHuman'], splitLine)

        if outProps[2] > 0.5:
            return
        lineDict['mHuman'] = outProps

        spec_ids = ['Pan', 'Pongo', 'Gorilla']

        for idS in spec_ids:
            outProps = checkFreq(specs[idS], splitLine)
            if outProps[2] > 0.5:
                return
            lineDict[idS] = outProps
    else:
        print(alt)
        return

    mHprops = lineDict['mHuman']
    mhAllele = mHprops.index(max(mHprops))

    alleles = set()
    for idS in spec_ids:
        specProps = lineDict[idS]
        alleles.add(specProps.index(max(specProps)))

    if len(alleles) > 1:
        return
    elif mhAllele in alleles:
        return
    else:
        isoLine = iso[(iso['chrom'] == isochr) & (iso['chrom_pos_1'] == int(pos))]
        if len(isoLine) == 0:
            print(chrnum, pos, lineDict, 'missing in isoform tbl')
            return
        site = list(isoLine['nuc_pos'])[0]
        refseq = list(isoLine['refcore'])[0]
        epLine = ep[(ep['refcore'] == refseq) & (ep['aa_pos'] == site)]
        if len(epLine) == 0:
            print(chrnum, pos, lineDict, 'missing in ep tbl')
            return

        strand = checkstrand(iso[iso['refcore'] == refseq])

        PosInCodon = site % 3
        if PosInCodon == 0:
            if strand == "+":
                codon = fa.fetch(fa.references[0], start=int(pos)-3, end=int(pos))
            else:
                codon = fa.fetch(fa.references[0], start=int(pos)-1, end=int(pos)+2)
        elif PosInCodon == 1:
            if strand == "+":
                codon = fa.fetch(fa.references[0], start=int(pos)-1, end=int(pos)+2)
            else:
                codon = fa.fetch(fa.references[0], start=int(pos)-3, end=int(pos))
        elif PosInCodon == 2:
            codon = fa.fetch(fa.references[0], start=int(pos)-2, end=int(pos)+1)
        else:
            print(pos, PosInCodon)

        if strand == '+':
            aa = seq.translate(codon)
            altCodon = list(codon)
            altCodon[PosInCodon] = alt
            altCodon = ''.join(altCodon)
            altAA = seq.translate(altCodon)
        else:
            aa = seq.translate(seq.reverse_complement(codon))
            altCodon = list(codon)
            altCodon[PosInCodon] = alt
            altCodon = ''.join(altCodon)
            altAA = seq.translate(seq.reverse_complement(altCodon))
        refEP = list(epLine[aa.lower()])[0]
        altEP = list(epLine[altAA.lower()])[0]
        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(pos), str(lineDict), codon, aa, str(refEP), altCodon, altAA, str(altEP)))


def main(vcffn, loc, chrnum):
    o = open('output.txt', 'w')
    ep = pd.read_csv('/mnt/d/adaptive_polymorphisms/snp_mapping/ep_data-46spp_hg19.txt', sep = "\t")
    iso = pd.read_csv('/mnt/d/adaptive_polymorphisms/snp_mapping/refseq2chrom_map_hg19.primary_isoforms.txt', sep = '\t')
    fa = pysam.FastaFile('/mnt/z/hg19_reference/chr{}.fa.gz'.format(chrnum))
    if vcffn.endswith('.gz'):
        readGZIPvcf(vcffn, loc, o, ep, iso, fa)
    else:
        readVCF(vcffn, loc, o, ep, iso, fa)
    o.close()

if __name__ == '__main__':
    #vcf_file = sys.argv[1]
    vcf_file = 'chr22_cds_primates.vcf.gz'
    loc = '../ep/'
    chrnum = '22'
    main(vcf_file, loc, chrnum)
