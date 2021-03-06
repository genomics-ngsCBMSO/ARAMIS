#!/bin/bash

ARAMIS_PATH='/usr/local/bin/'

usage()
{
    echo "Usage: $0 -a <ASSEMBLY> -b1 <ILUALIGNMENT> -b2 <PBALIGNMENT> -t <TARGETS> -p <PREFIX>"
    echo "  -a | --assembly_file             PacBio assembly in Fasta format"
    echo " -b1 | --illumina_alignment_file   input illumina alignment in BAM (sorted) format"
    echo " -b2 | --pacbio_alignment_file     input pacbio alignment in BAM (sorted) format"
    echo "  -p | --prefix                    Prefix for output files"
    echo "  -t | --targets_file              Target file generated by correction.sh"
    echo "  -h | --help                      Display help"

}

ASSEMBLY=''
ILUALIGNMENT=''
PBALIGNMENT=''
TARGETS=''
PREFIX='standard'


# display usage if 
	if [ $# -lt 6 ]; then
	    echo "You must provide at least the assembly and alignment files"
	    usage
	    exit 2
	fi


while [[ "$1" > 0 ]]
do
  case $1 in
   -a| --assembly_file)
    shift
    ASSEMBLY=$1
    shift
  ;;
   -b1| --illumina_alignment_file)
    shift
    ILUALIGNMENT=$1
    shift
  ;;
   -b2| --pacbio_alignment_file)
    shift
    PBALIGNMENT=$1
    shift
  ;;
   -t| --targets_file)
    shift
    TARGETS=$1
    shift
  ;;
    -p|--prefix )
    shift
    PREFIX=$1
    shift
  ;; 
    -h|--help)
    usage
    exit
  ;;
  *)
    echo "Wrong parameter!!!"
  ;;
  esac
done


echo -e "\n"

echo "- The file $ASSEMBLY will be used as the reference assembly.

- The $ILUALIGNMENT will be used as Illumina alignment file.

- The $PBALIGNMENT will be used as Pacbio alignment file.

- You selected the file $TARGETS generated by correction.sh"


echo -e "\n"
echo "Let's start calculating coverage !! (This could take a while)"

mkdir indel_information &>/dev/null
cd indel_information/

ln -s ../$ASSEMBLY &>/dev/null
ln -s ../$ILUALIGNMENT &>/dev/null
ln -s ../$PBALIGNMENT &>/dev/null
ln -s ../$TARGETS &>/dev/null

##LET'S START CALCULATING COVERAGE

samtools index $ILUALIGNMENT
samtools index $PBALIGNMENT 

igvtools count -w 1 --bases --minMapQuality 1 $ILUALIGNMENT ${PREFIX}_Illumina_coverage.wig $ASSEMBLY &>/dev/null

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE IGVTOOLS: Illumina coverage!!"
		exit 1
	fi 


igvtools count -w 1 --bases --minMapQuality 1 $PBALIGNMENT ${PREFIX}_PacBio_coverage.wig $ASSEMBLY &>/dev/null

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE IGVTOOLS: PacBio coverage!!"
		exit 1
	fi 

echo -e "\n"
echo "Coverage calculation DONE !!"

#NOW DETECT ALL HOMOPOLYMERS IN THE ASSEMBLY 

#Motifs.fa: Fasta file with the motifs to search: homopolymers with at least 3 G/C/T/A

echo -e "\n"
echo "Finding homopolymers tracks !!"

echo -e '>Cs\nC{2,}\n\n>Gs\nG{2,}\n\n>As\nA{2,}\n\n>Ts\nT{2,}\n' > motifs.fa

seqkit locate -P -f motifs.fa $ASSEMBLY | sort -n -k 5 | awk 'BEGIN{FS=OFS="\t"}{print $1,$5,$6,$4,$2,$7,$6-$5+1}' > ${PREFIX}_homopolymers.txt


echo -e "\n"
echo "All homopolymers detected !!"


#CONTINUE CALCULATING GC SKEW

echo -e "\n"
echo "Let's calculate GC-SKEW !!"

gc_skew.py -f $ASSEMBLY -s 1000 -w 1000 --no-plot > ${PREFIX}_gc_skew.txt  

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE GC SKEW CALCULATOR !!"
		exit 1
	fi 

echo -e "\n"
echo "GC-SKEW calculated !!"


#COMBINING ALL INFORMATION

echo -e "\n"
echo "Let's combine all information and generate some plots !!"

echo -e "Indels\t$TARGETS\nIllumina_coverage\t${PREFIX}_Illumina_coverage.wig\nPacBio_coverage\t${PREFIX}_PacBio_coverage.wig\nHomopolymers_file\t${PREFIX}_homopolymers.txt\nGCSkew_file\t${PREFIX}_gc_skew.txt" > list

samtools faidx $ASSEMBLY
cut -f1-2 ${ASSEMBLY}.fai > list_chro

python3 ${ARAMIS_PATH}combine_info.py list ${PREFIX}_allinfo.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE A PYTHON3 SCRIPT !!"
		exit 1
	fi 

Rscript --vanilla ${ARAMIS_PATH}plot_generation.R -f ${PREFIX}_allinfo.txt -l list_chro -c tempcoverage.txt -g ${PREFIX}_gc_skew.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE A R SCRIPT !!"
		exit 1
	fi 

rm tempcoverage.txt list_chro list motifs.fa
find -type l -delete

echo -e "\n"
echo "ALL DONE. Now check the results !!"


echo "                            /   |
 _                        )      ((   ))     (
(@)                      /||      ))_((     /||
|-|                     / | |    (/||/|)   / | |                      (@)
| | -------------------/--|-voV---|'|'/--Vov-|--|---------------------|-|
|-|                         '^'   (o o)  '^'                          | |
| |                               '|Y/'                               |-|
|-|                                                                   | |
| |                         That's all folks!!!                       |-|
|-|                                                                   | |
| |  Eva Sacristán Horcajada                                          |-|
| |            Sandra González de la Fuente                           | |
|-|                         Ramón Peiró Pastor                        |-|
|_|___________________________________________________________________| |
(@)              l   /| /         ( (       | /|   l                '||-|
                 l /   V           | |       V   | l                  (@)
                 l/                _) )_          |I
                                   '| /'
"






