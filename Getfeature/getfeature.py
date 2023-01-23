
import re, os, sys
import re
from collections import Counter
import Bio.SeqIO as Seq  #导入Bio包
from feamodule import fickett
import argparse as agp
from nucleotide_feature1 import get_RNA_feature
from protein_feature3 import get_protein_feature
from nucleotide_feature2 import get_fhCTD_feature


def get_feature():

    parser = agp.ArgumentParser()
    parser.add_argument('-i', '--nucleotide_file', help="the input FASTA file of nucleotide sequence")
    #parser.add_argument('-ip', '--protein_file', help="the input FASTA file of protein sequence")
    parser.add_argument('-o', '--outfile', help="output file")
    args = parser.parse_args()

    seq_ntfile=args.nucleotide_file
    outfile=args.outfile

    f = open(outfile, 'w')
    with open(seq_ntfile, 'r') as r:
        lines = r.readlines()
    for seq in lines:
        feature_nt=get_RNA_feature(seq)
        feature_hfCTD=get_fhCTD_feature(seq)

        lister = []
        seq_new = ''
        for j in range(0, len(seq) - 2, 3):
            if (seq[j] == 'T' and seq[j + 1] == 'T' and seq[j + 2] == 'T') or (seq[j] == 'T' and seq[j + 1] == 'T' and seq[j + 2] == 'C'):
                seq_new = seq_new + 'F'
            elif (seq[j] == 'T' and seq[j + 1] == 'T' and seq[j + 2] == 'A') or (seq[j] == 'T' and seq[j + 1] == 'T' and seq[j + 2] == 'G') or (seq[j] == 'C' and seq[j + 1] == 'T'):
                seq_new = seq_new + 'L'
            elif (seq[j] == 'A' and seq[j + 1] == 'T' and seq[j + 2] == 'T') or (seq[j] == 'A' and seq[j + 1] == 'T' and seq[j + 2] == 'C') or (seq[j] == 'A' and seq[j + 1] == 'T' and seq[j + 2] == 'A'):
                seq_new = seq_new + 'I'
            elif (seq[j] == 'A' and seq[j + 1] == 'T' and seq[j + 2] == 'G'):
                seq_new = seq_new + 'M'
            elif (seq[j] == 'G' and seq[j + 1] == 'T'):
                seq_new = seq_new + 'V'
            elif (seq[j] == 'T' and seq[j + 1] == 'C') or (seq[j] == 'A' and seq[j + 1] == 'G' and seq[j + 2] == 'T') or (seq[j] == 'A' and seq[j + 1] == 'G' and seq[j + 2] == 'C'):
                seq_new = seq_new + 'S'
            elif (seq[j] == 'C' and seq[j + 1] == 'C'):
                seq_new = seq_new + 'P'
            elif (seq[j] == 'A' and seq[j + 1] == 'C'):
                seq_new = seq_new + 'T'
            elif (seq[j] == 'G' and seq[j + 1] == 'C'):
                seq_new = seq_new + 'A'
            elif (seq[j] == 'T' and seq[j + 1] == 'A' and seq[j + 2] == 'T') or (seq[j] == 'T' and seq[j + 1] == 'A' and seq[j + 2] == 'C'):
                seq_new = seq_new + 'Y'
            elif (seq[j] == 'C' and seq[j + 1] == 'A' and seq[j + 2] == 'T') or (seq[j] == 'C' and seq[j + 1] == 'A' and seq[j + 2] == 'C'):
                seq_new = seq_new + 'H'
            elif (seq[j] == 'C' and seq[j + 1] == 'A' and seq[j + 2] == 'A') or (seq[j] == 'C' and seq[j + 1] == 'A' and seq[j + 2] == 'G'):
                seq_new = seq_new + 'Q'
            elif (seq[j] == 'A' and seq[j + 1] == 'A' and seq[j + 2] == 'T') or (seq[j] == 'A' and seq[j + 1] == 'A' and seq[j + 2] == 'C'):
                seq_new = seq_new + 'N'
            elif (seq[j] == 'A' and seq[j + 1] == 'A' and seq[j + 2] == 'A') or (seq[j] == 'A' and seq[j + 1] == 'A' and seq[j + 2] == 'G'):
                seq_new = seq_new + 'K'
            elif (seq[j] == 'G' and seq[j + 1] == 'A' and seq[j + 2] == 'T') or (seq[j] == 'G' and seq[j + 1] == 'A' and seq[j + 2] == 'C'):
                seq_new = seq_new  + 'D'
            elif (seq[j] == 'G' and seq[j + 1] == 'A' and seq[j + 2] == 'A') or (seq[j] == 'G' and seq[j + 1] == 'A' and seq[j + 2] == 'G'):
                seq_new  = seq_new  + 'E'
            elif (seq[j] == 'T' and seq[j + 1] == 'G' and seq[j + 2] == 'T') or (seq[j] == 'T' and seq[j + 1] == 'G' and seq[j + 2] == 'C'):
                seq_new  = seq_new  + 'C'
            elif (seq[j] == 'T' and seq[j + 1] == 'G' and seq[j + 2] == 'G'):
                seq_new  = seq_new  + 'W'
            elif (seq[j] == 'C' and seq[j + 1] == 'G') or (seq[j] == 'A' and seq[j + 1] == 'G' and seq[j + 2] == 'A') or (seq[j] == 'A' and seq[j + 1] == 'G' and seq[j + 2] == 'G'):
                seq_new  = seq_new  + 'R'
            elif (seq[j] == 'G' and seq[j + 1] == 'G'):
                seq_new  = seq_new  + 'G'

        lister.append(seq_new )
        print(lister)
        seq_aa="".join(map(str, lister))
        print(seq_aa)

        feature_aa=get_protein_feature(seq_aa)


        feature = feature_nt + feature_aa + feature_hfCTD
        s = str(feature).replace('[', '').replace(']', '')
        f.write(s+'\n')
get_feature()

