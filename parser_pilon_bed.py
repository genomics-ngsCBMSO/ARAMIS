#!/usr/bin/env python3
# Imports

from __future__ import division
from collections import defaultdict
import sys, os
import math
import re

syntax = '''----------------------------------------------------------------------------
	Description: Script to parse Pilon changes output file to a bed file.

        #Usage: python3 parser_pilon_bed.py <pilon.changes_file>

-----------------------------------------------------------------------------'''

if len(sys.argv) != 2:
  print (syntax)
  sys.exit()

pil_file = sys.argv[1]
out_file = pil_file + ".bed"

inhandle = open(pil_file, 'r')
outhandle = open(out_file, 'w')

print ("    Opening and reading ", pil_file, "as Pilon changes file ....")
lines = inhandle.readlines()
for line in lines:
  newline = re.sub('(:)', ' ', line)
  #print(newline)
  newline = newline.rstrip("\n")
  chunk = newline.split(" ")
  chro = chunk[0]
  #ori_pos = chunk[1]
  if "-" in chunk[1]:
    ori_pos = chunk[1].split("-")
    ori_pos_1 = ori_pos[0]
    ori_pos_2 = ori_pos[1]
  else: 
    ori_pos_1 = chunk[1]
    ori_pos_2 = "NA"
  ori_base = chunk[4]
  fin_base = chunk[5]
  if "-" in chunk[3]:
    fin_pos = chunk[3].split("-")
    fin_pos_1 = fin_pos[0]
    fin_pos_2 = fin_pos[1]
  else:
    fin_pos_1 = chunk[3]
    fin_pos_2 = "NA"
  data = str(chro) + "\t" + str(ori_pos_1) + "\t" + str(ori_pos_2) + "\t" + str(fin_pos_1) + "\t" + str(fin_pos_2) + "\t" + str(ori_base) + "\t" + str(fin_base)  + "\n"
  outhandle.write(data)
    
inhandle.close()
outhandle.close()

print("    Output file generated: ", out_file)
print("ALL DONE!!!")

