#!/usr/bin/env python
# _*_coding:utf-8_*_
import os,sys
from feamodule import CTD
import numpy
import Bio.SeqIO as Seq
from feamodule import fickett
from feamodule import FrameKmer


def coding_nocoding_potential(input_file):
	coding={}
	noncoding={}
	for line in open(input_file).readlines():
		fields = line.split()
		if fields[0] == 'hexamer':continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] = float(fields[2])
	return coding,noncoding


def get_fhCTD_feature(seq):
	hex_file='ath2300_hexamer_tab.tsv'  #hexamer 的频率表要根据数据集的不同而更改
	coding,noncoding = coding_nocoding_potential(hex_file)
	fickett_fe = fickett.fickett_value(seq)
	hexamer = FrameKmer.kmer_ratio(seq,6,3,coding,noncoding)
	A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = CTD.CTD(seq)

	feature=[fickett_fe,hexamer,float(A),float(T),float(G),float(C),float(AT),float(AG),float(AC),float(TG),float(TC),float(GC),
			 float(A0),float(A1),float(A2),float(A3),float(A4),float(T0),float(T1),float(T2),float(T3),float(T4),
			 float(G0),float(G1),float(G2),float(G3),float(G4),float(C0),float(C1),float(C2),float(C3),float(C4)]
	return feature

