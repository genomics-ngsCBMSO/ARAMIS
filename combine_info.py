#!/usr/bin/env python3
# Imports

from __future__ import division
from collections import defaultdict
import sys, os
import math
import csv
import time
import fileinput
import subprocess

syntax = '''----------------------------------------------------------------------------
	Description: Script to combine all the information about the corrected errors 
	obtained in previous steps from the script indel_analysis.sh

        #Usage: python3 combine_info.py <file_list> <out_file>

        file_list : Tabular file indicating all the results files generated in previous steps
	The list MUST have the following structure:

		Indels			<indels_file>
		Illumina_coverage	<Illumina_coverage_file>
		PacBio_coverage		<PacBio_coverage_file>
		Homopolymers_file	<Homopolymers_file>
		GCSkew_file		<gc_skew_file>

-----------------------------------------------------------------------------------------'''

if len(sys.argv) != 3:
  print (syntax)
  sys.exit()

for line in fileinput.input(sys.argv[1]):
  fileline = line.strip()
  filename = fileline.split("\t")
  #print(filename)
  if filename[0] == "Indels":  
    pos_file = filename[1]
    inhandle1 = open(pos_file, 'r')
  elif filename[0] == "Illumina_coverage":  
    covilu_file = filename[1]
    inhandle2 = open(covilu_file, 'r')
  elif filename[0] == "PacBio_coverage":
    covpb_file = filename[1]
    inhandle3 = open(covpb_file, 'r')
  elif filename[0] == "Homopolymers_file":
    hom_file = filename[1]
    inhandle4 = open(hom_file, 'r')
  elif filename[0] == "GCSkew_file":
    gcs_file = filename[1]
  else:
    print("\t ERROR: The list of files has an incorrect format!!")
    print(syntax)
    sys.exit()
   
out_file1 = sys.argv[2]
#out_file2 = sys.argv[3]

outhandle1 = open(out_file1, 'w')
outhandle2 = open("tempcoverage.txt", 'w')
outhandle3 = open("stats.txt", 'w')

pos_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
covilu_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
covpb_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
hom_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
pal_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))


#Read Indel Positions file

print("\n Reading Indel positions \n")
indel_count = 0
lines = inhandle1.readlines()
for line in lines:
  if "#" in line:
    continue
  else: 
    line = line.rstrip("\n")
    chunk = line.split("\t")
    chro = chunk[0]
    pos = int(chunk[1])
    change = chunk[8]
    pos_dic[chro][pos]['change'] = change
    indel_count = indel_count + 1     

inhandle1.close() 

print("\t","Number of indels = ", indel_count, "\n")

#Read Illumina IgvTools count output file

print(" Reading Illumina coverage information \n")
lines = inhandle2.readlines()
for line in lines:
  line = line.rstrip("\n")
  if line[0] == "#" or "track type" in line:
    continue
  elif "variableStep" in line:
    string = line.split(" ")
    substring = string[1].split("=")
    chro = substring[1]
    print ("\t", "Parsing", chro, "\n")
  else: 
    chunk = line.split("\t") 
    pos = int(chunk[0])
    cov  = float(chunk[1]) + float(chunk[2]) + float(chunk[3]) + float(chunk[4])
    dels = float(chunk[6])
    ins = float(chunk[7])
    covilu_dic[chro][pos]['cov'] = cov
    covilu_dic[chro][pos]['ins'] = ins
    covilu_dic[chro][pos]['dels'] = dels

inhandle2.close()

#Read PacBio IgvTools count output file

print(" Reading PacBio coverage information \n")
lines = inhandle3.readlines()
for line in lines:
  line = line.rstrip("\n")
  if line[0] == "#" or "track type" in line:
    continue
  elif "variableStep" in line:
    string = line.split(" ")
    substring = string[1].split("=")
    chro = substring[1]
    print ("\t", "Parsing", chro, "\n")
  else: 
    chunk = line.split("\t") 
    pos = int(chunk[0])
    cov  = float(chunk[1]) + float(chunk[2]) + float(chunk[3]) + float(chunk[4])
    dels = float(chunk[6])
    ins = float(chunk[7])
    covpb_dic[chro][pos]['cov'] = cov
    covpb_dic[chro][pos]['ins'] = ins
    covpb_dic[chro][pos]['dels'] = dels

inhandle3.close()

#Read Homopolymers file

print(" Reading Hompolymers information \n")
hom_count,totalA,totalT,totalC,totalG = 0,0,0,0,0
lines = inhandle4.readlines()
for line in lines:
  line = line.rstrip("\n")
  if "seqID" in line:
    continue
  else:
    chunk = line.split("\t")
    chro = chunk[0]
    pos = int(chunk[1]) - 1
    end = int(chunk[2]) - 1
    base = chunk[4]
    lgth = int(chunk[6])
    if base == "As":
      base = "A"
      totalA = totalA + 1
    elif base == "Ts":
      base = "T"
      totalT = totalT + 1
    elif base == "Cs":
      base = "C"
      totalC = totalC + 1
    else:
      base = "G"
      totalG = totalG + 1
    hom_dic[chro][pos]['base'] = base
    hom_dic[chro][pos]['lgth'] = lgth
    hom_dic[chro][pos]['end'] = end
    hom_count = hom_count + 1

inhandle4.close()

print("\t", "Number of homopolymers: ", hom_count, "\n")

#Combine information

header1 = "Chromosome" + "\t" + "Position" + "\t" + "Change" + "\t" + "Illumina_cov" + "\t" + "Illumina_Indel_fraction" + "\t" + "PacBio_cov" + "\t" + "PacBio_Indel_fraction" +"\t" + "Homopolymer" + "\t" + "Type" + "\t" + "Length" + "\n"
outhandle1.write(header1)

header2 = "Chromosome" + "\t" + "Position" + "\t" + "Coverage" + "\t" + "Indel_fraction" + "\t" + "Platform" + "\n"
outhandle2.write(header2)

print(" Merging all information and calculating basic stats \n")
print(" This could take a while... There are a lot of lines!\n")
homyes_count, A_count, G_count, C_count, T_count = 0,0,0,0,0
for chro in sorted(pos_dic.keys()):
  for pos in sorted(pos_dic[chro].keys()):
    ilucov = covilu_dic[chro][pos]['cov']
    pbcov = covpb_dic[chro][pos]['cov']
    change = pos_dic[chro][pos]['change']
    if '-' in change:
      pos = pos + 1
      iludels = covilu_dic[chro][pos]['dels']
      ilufrac = float(iludels) / float(ilucov)
      pbdels = covpb_dic[chro][pos]['dels']
      pbfrac = float(pbdels) / float(pbcov)
      pos = pos - 1
    else:
      iluins = covilu_dic[chro][pos]['ins']
      ilufrac = float(iluins) / float(ilucov)
      pbins = covpb_dic[chro][pos]['ins']
      pbfrac = float(pbins) / float(pbcov)
    if ilufrac > 1 or pbfrac > 1:
      ilufrac,pbfrac = 1,1
    data_covilu = str(chro) + "\t" + str(pos) + "\t" + str(ilucov) + "\t" + str(ilufrac) + "\t" + "Illumina" + "\n"
    outhandle2.write(data_covilu)
    data_covpb = str(chro) + "\t" + str(pos) + "\t" + str(pbcov) + "\t" + str(pbfrac) + "\t" + "PacBio" + "\n"
    outhandle2.write(data_covpb)
    if pos in sorted(hom_dic[chro].keys()):
      hom = "YES"
      homyes_count = homyes_count + 1
      base = hom_dic[chro][pos]['base']
      if base == "A":
        A_count = A_count + 1
      elif base == "T":
        T_count = T_count + 1
      elif base == "C":
        C_count = C_count + 1
      else:
        G_count = G_count + 1
      lgth = hom_dic[chro][pos]['lgth']
    else:
      hom = "NO"
      base = "NA"
      lgth = "NA"
    data = str(chro) + "\t" + str(pos) + "\t" + str(change) + "\t" + str(ilucov) + "\t" + str(ilufrac) + "\t" + str(pbcov) + "\t" + str(pbfrac) + "\t" + str(hom) + "\t" + str(base) + "\t" + str(lgth) + "\n" 
    outhandle1.write(data)

#Generate Stats File

col1 = ["Total Homopolymers", "Total Indels", "Indels in Homopolymers", "% of Homopolymers Affected","A homopolymers with indel", "% of A Homopolymers Affected", "T homopolymers with indel", "% of T Homopolymers Affected","G homopolymers with indel", "% of G Homopolymers Affected","C homopolymers with indel", "% of C Homopolymers Affected"]
col2 = [hom_count, indel_count, homyes_count, ((homyes_count/hom_count)*100), A_count, ((A_count/totalA)*100), T_count, ((T_count/totalT)*100), G_count, ((G_count/totalA)*100), C_count, ((C_count/totalA)*100)]

writer = csv.writer(outhandle3, delimiter='\t')
writer.writerows(zip(col1,col2))

outhandle1.close()
outhandle2.close()
outhandle3.close()

#That's all

print(" Output files generated:", out_file1,"+ stats.txt + plots")
