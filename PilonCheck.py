#!/usr/bin/env python3
# Imports
from __future__ import division
from collections import defaultdict
import sys, os
import math

syntax = '''----------------------------------------------------------------------------
	Descripcion: El script genera un archivo en formato PacBio utilities con las posiciones bad de PacBio utilities que se han generado por presencia
	de mas de una opcion de correccion. Las posiciones correctas a cambiar se extraen del changes de Pilon(que previamente se ha pasado a formato bed con parser_pilon_bed.pl)
	Ademas genera un archivo con aquellos posibles cambios de PacBio utilities que no se encuentra en Pilon).
	Usage: python PilonCheck.py <Pilon.changes.bed> <targets_onlybad.txt> <pilon_common.txt> <pilon_not_common.txt> <not_warning_file.txt> <warning_file.txt>
-----------------------------------------------------------------------------'''
if len(sys.argv) != 7:
        print (syntax)
        sys.exit()
        
#--Parameters
changesPilon_1 = sys.argv[1] 
targetsPacBioutilities_1 = sys.argv[2]
output_file = sys.argv[3]
output2_file = sys.argv[4]
output3_file = sys.argv[5]
output4_file = sys.argv[6]

inhandle1 = open(changesPilon_1 , 'r')
inhandle2 = open(targetsPacBioutilities_1 , 'r')
outhandle1 = open(output_file , 'w')
outhandle2 = open(output2_file , 'w')
outhandle3 = open(output3_file , 'w')
outhandle4 = open(output4_file , 'w')

# lee el archivo de Pilon. 
pilon_dic = defaultdict(lambda: defaultdict(lambda: defaultdict()))
#header = 
#outhandle.write(header)

lines1 = inhandle1.readlines()
for line1 in lines1:
  if "#" in line1:
    continue
  else:
    line1 = line1.rstrip("\n")
    chunks1 = line1.split("\t")
    name1 = chunks1[0]
    inicio1 = int(chunks1[1])-1
    pos1 = chunks1[5]
    pos2 = chunks1[6]
    lengthBasePos1 = len(pos1)
    lengthBasePos2 = len(pos2)
    pilon_dic[name1][inicio1]['line1'] = line1
    pilon_dic[name1][inicio1]['pos1'] = pos1
    pilon_dic[name1][inicio1]['pos2'] = pos2
    pilon_dic[name1][inicio1]['lengthBasePos1'] = lengthBasePos1
    pilon_dic[name1][inicio1]['lengthBasePos2'] = lengthBasePos2
inhandle1.close()


# lee el archivo de PacBio utilities. 

lines2 = inhandle2.readlines()
for line2 in lines2:
  if "#" in line2:
    outhandle1.write(line2)
  else:
    line2 = line2.rstrip("\n")
    chunks2 = line2.split("\t")
    name2 = chunks2[0]
    inicio2 = int(chunks2[1])
    haplotipo = chunks2[8]
    if "," in haplotipo:
      haplotipo = haplotipo
    try:
      if pilon_dic[name2][inicio2]['pos1'] == '.':
        if haplotipo:
           outhandle4.write(line2 + "\t" + "Warning, you should check this target!!" + "\n")
        else:
           outhandle3.write(line2 + "\t" + "Non warning, good target!!" + "\n")
        longitudCambio = str(pilon_dic[name2][inicio2]['lengthBasePos2'])
        cambio = ("+" + str(pilon_dic[name2][inicio2]['lengthBasePos2']) + str(pilon_dic[name2][inicio2]['pos2']))       
        outhandle1.write(str(name2) + "\t" + str(inicio2) + "\t" + str(chunks2[2]) + "\t" +  str(chunks2[3]) + "\t" + \
        str(chunks2[4]) + "\t" +  str(chunks2[5]) + "\t" +  str(chunks2[6]) + "\t" + str(longitudCambio) + "\t" + \
        str(cambio) + "\t" + str(chunks2[9]) + "\t" +  str(chunks2[10]) + "\n")
      else:
        if haplotipo:
           outhandle4.write(line2 + "\t" + "Warning, you should check this target!!" + "\n")
        else:
           outhandle3.write(line2 + "\t" + "Non warning, good target!!" + "\n")
        longitudCambio = str(pilon_dic[name2][inicio2]['lengthBasePos1'])
        cambio = ("-" + str(pilon_dic[name2][inicio2]['lengthBasePos1']) + str(pilon_dic[name2][inicio2]['pos1']))
        outhandle1.write(str(name2) + "\t" + str(inicio2) + "\t" + str(chunks2[2]) + "\t" +  str(chunks2[3]) + "\t" + \
        str(chunks2[4]) + "\t" +  str(chunks2[5]) + "\t" +  str(chunks2[6]) + "\t" + '-' + str(longitudCambio) + "\t" +\
        str(cambio) + "\t" + str(chunks2[9]) + "\t" +  str(chunks2[10]) + "\n")
    except KeyError:
      outhandle2.write(line2 + "\n")

inhandle2.close()

print("ALL DONE!!!")