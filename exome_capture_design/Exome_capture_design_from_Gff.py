#! /usr/bin/env python3

import numpy as np
import os
import sys
import re
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

Usage = """
Gff_editscript.py, v. 1.0, last updated 12 April 2019, Phred M. Benham UC Berkeley.
Gff_editscript.py parses GFF file to produce fasta file of exons (or other regions) passing user-specified filters.
Can specify different regions of genomic intervals in gff, interval length, acceptable gc content, add flanking regions to intervals, etc.
Input: python3 Gff_editscript.py -g <gff file for reference> -f <fasta/fna/fa file for reference genome>.
Output: Fasta file of sequences for all exons passing filtering steps.
Requires: bedtools (https://bedtools.readthedocs.io/en/latest/) and RepeatMasker (http://www.repeatmasker.org).
"""

parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-g","--gff", required=True, help="input gff files")
parser.add_argument("-f","--fasta", required=True, help="input fna/fasta files")
parser.add_argument("-l","--genome", required=True, help="input genome file (contig lengths)")
args = parser.parse_args()

#print(args.gff)
#print(args.fasta)

#input options to be user specified
MyGFF = args.gff
MyFa = args.fasta
MyGenFile = args.genome
region = 'exon'
region_length = 150
minGC = 0.3
maxGC = 0.7
flank = 15



#Optional step to remove loci that are uncharacterized in reference gff.
#In the case of Zonotrichia albicollis genome these were identified as expressed transcripts, but are unannotated.
print("remove uncharacterized loci")
os.system("grep -v product=uncharacterized %s >unchar.removed.gff" % MyGFF)


newGFF = "unchar.removed.gff"
#use bedtools to annotate original .gff file with info on gcc content, feature length, etc.
print("annontate gff with GC content")
os.system("bedtools nuc -fi %s -bed %s >WTSP_GCannontate.gff" % (MyFa, newGFF))

#Functions to return new GFF restricted to
#exons greater than region_length and with GC content between minGC-maxGC

print("filtering gff based on exon, exon length, GC content")
Gff_input = open("WTSP_GCannontate.gff", 'r')
OutFile = open("WTSP_finalexons.gff", 'w')

#create dictionary with scaffold name_start_end position to only include single transcript variant and not redundant sequence info.
Unique_exon = {}
for Line in Gff_input:
	#skip first line
	if Line[0] == '#':
		OutFile.write(Line)
	else:
		Line=Line.strip('\n')
		ColumnList = Line.split('\t')
		
		#only print line for exons greater than 100bp length and with GC content between 0.3-0.7%
		if ColumnList[2] == region and int(ColumnList[17]) > region_length and float(ColumnList[10])<=maxGC and float(ColumnList[10])>=minGC:
			SeqKey = ColumnList[0] + '_' + ColumnList[3] + '_' + ColumnList[4]
			Unique_exon[SeqKey] =''
			Unique_exon[SeqKey] += Line

#loop through dictionary to print lines of retained exons to new gff file.
#also calculate the total number of exons retained and sum length in bp of these regions.			

print("writing outfile with exons extended by 15bp")

Genome_input = open(MyGenFile,'r')
Genome_dict = {}
for Line in Genome_input:
	Line = Line.strip('\n')
	LineList = Line.split('\t')
	ChrKey = LineList[0]
	Genome_dict[ChrKey] = ''
	Genome_dict[ChrKey] += LineList[1]
	
ExonNumber = 0
TotalSequenceLength = 0
for key in Unique_exon:
	ColumnList2=Unique_exon[key].split('\t')
	startbp = str(int(ColumnList2[3])-flank)
	endbp = str(int(ColumnList2[4])+flank)
	
	beginlist = '\t'.join(ColumnList2[0:3])
	endlist = '\t'.join(ColumnList2[5:18])
	
	if int(startbp) <= 0 and int(Genome_dict[ColumnList2[0]]) < int(endbp):
		OutFile.write(beginlist + '\t' + '1' + '\t' + Genome_dict[ColumnList2[0]] + '\t' + endlist + '\n')
		
	elif int(startbp) <= 0:
		#print line with column 3 = 1
		OutFile.write(beginlist + '\t' + '1' + '\t' + endbp + '\t' + endlist + '\n')
	
	elif int(Genome_dict[ColumnList2[0]]) < int(endbp):
		OutFile.write(beginlist + '\t' + startbp + '\t' + Genome_dict[ColumnList2[0]] + '\t' + endlist + '\n')

	else:
		#print line normal	
		OutFile.write(beginlist + '\t' + startbp + '\t' + endbp+ '\t' + endlist + '\n')
	
	ExonNumber += 1
	TotalSequenceLength += int(ColumnList2[17])

temp = sys.stdout
sys.stdout = open('Gff_edit.log.txt', 'w')
print("Exons information following GC, length filters:")
print("Number of exons retained:\t%s" % str(ExonNumber))
print("Sequence length of retained exons:\t%s\n" % str(TotalSequenceLength))


OutFile.close()
Gff_input.close()	

FinGff = open("WTSP_finalexons.gff", 'r')
FaHead_dict = {}
for Line in FinGff:
	#skip first line
	Line=Line.strip('\n')
	if Line[0] == '#':
		print('oops')
	else:	
		search_term ='^(.+)\t\w.+\texon\t(\d+)\t(\d+)\t.+(gene=.+?)\;.+'
		Results = re.search(search_term, Line)
		print(Results)
		#NewVal = int(Results.group(2)) - 1
		#Fa_head = '>' + Results.group(1) + ':' + str(NewVal) #+ '-' + Results.group(3) 
		#GeneID = '-' + Results.group(3) + ':' + Results.group(4)
		#FaHead_dict[Fa_head]=''
		#FaHead_dict[Fa_head] += GeneID

FinGff.close()	


sys.stderr.write("writing fasta file")
os.system("sort -k1,1 -k4,4n WTSP_finalexons.gff >WTSP_finalexons.sort.gff")
os.system("bedtools merge -i WTSP_finalexons.sort.gff >WTSP_finalexons.merge.gff")
os.system("bedtools getfasta -fi %s -bed WTSP_finalexons.merge.gff -fo WTSP_finalexons.fa" % MyFa) 


#Run RepeatMasker on resulting fasta file
sys.stderr.write("Running repeat masker...")
RMcommand = "repeatmasker -s -species chicken  WTSP_finalexons.fa" 
os.system(RMcommand)

#really important awk command that I forgot the point of?
awkCommand = """awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\\n"$0} else {printf $0}}}'  WTSP_finalexons.fa.masked > WTSP_finalexons.masked.fa"""
os.system(awkCommand)

#final edit of fasta file to remove any sequences with repeats, etc. masked by NNNs 
#or sequences shorter than 1000bp from final file.
sys.stderr.write("final sequence editing...")
Infile = "WTSP_finalexons.masked.fa"
Outfile = "WTSP_exons_final.fa" 

Fa_input = open(Infile, 'r')
OutFile = open(Outfile, 'w')

Exon_length = []

#Last set of functions to create final fasta file of exons
# also calculate number of exons and summary stats from exons. 
ref_dict={}
for Line in Fa_input:
	Line = Line.strip('\n')
	if '>' in Line:
		Line2 = Line.split('-')
		seq = Line2[0]
		ref_dict[seq]=''
		
	else:
		ref_dict[seq] += Line

seq_remove = 0
seq_keep = 0
seq_length = 0
for seq in ref_dict:
	UnmaskedString = re.findall(r'[ACTG]+', ref_dict[seq])
	LongSeq = max(UnmaskedString)
	if len(LongSeq) >= 180:
		OutFile.write(seq + FaHead_dict[seq] + '\n' + LongSeq + '\n')
		seq_keep += 1
		seq_length+=len(LongSeq)
		Exon_length.append(len(LongSeq))
		
print("Final exon length and numbers:")		
print("Number of retained exons:\t%s" % str(seq_keep))
print("Total length of final exons:\t%s" % str(seq_length))

OutFile.close()
Fa_input.close()

#exon summary stats and histogram

Exon_length = np.asarray(Exon_length)
exon_length_mean = np.mean(Exon_length)
exon_length_std = np.nanstd(Exon_length)
exon_max = np.nanmax(Exon_length)
print("Mean Length:", exon_length_mean, "Sd:", exon_length_std, "Max exon length:", exon_max, sep='\t')
sys.stdout.close()
sys.stdout = temp


fig, ax = plt.subplots(figsize=(9, 7))
sns.despine(ax=ax, offset=10)
ax.hist(Exon_length, bins=50, edgecolor='black')
ax.set_xlabel('Exon length')
ax.set_ylabel('Number of exons')
ax.set_title('Distribution of exon lengths')

plt.savefig('ExonLengthHistogram.png')