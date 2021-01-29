#!/bin/bash

ARAMIS_PATH='~/'
PICARD_PATH='~/picard.jar'
PILON_PATH='~/pilon-1.22.jar'

usage()
{
    echo "Usage: $0 -a <ASSEMBLY_FILE> -b <ALIGNMENT_FILE> -i <INDEL_FRAC> -p <PREFIX>"
    echo -e "\n"
    echo "In case you want to perform manual correction of warning targets, please run Usage2"
    echo "Usage2: $0 -a <ASSEMBLY_FILE> -w  <WARNINGS> -p <PREFIX>"
    echo -e "\n"
    echo "  -a | --assembly_file    PacBio assembly in Fasta format (REQUIRED)"
    echo "  -b | --alignment_file   Input alignment in BAM (sorted) format (REQUIRED)"
    echo "  -i | --indel_fraction   Do not report indels at a position if the fraction of reads containing them is below FLOAT"
    echo " -pt | --picard_path      If picard.jar is not in the path /home/$USER, please indicate the correct path with this option"  
    echo " -pl | --pilon_path       If pilon.jar is not in the path /home/$USER, please indicate the correct path with this option" 
    echo "  -w | --warnings         Alternative pipeline after manual correction of indels flags as warning."
    echo "  -p | --prefix           Prefix for output files"
    echo "  -h | --help             Display help"

}

ASSEMBLY=''
ALIGNMENT=''
INDEL_FRAC=0.5
WARNINGS=''
PREFIX='standard'

# display usage if 
	if [ $# -lt 4 ]; then
	    echo "You must provided at least the assembly and alignment files"
	    usage
	    exit 2
	fi

#Get parameters

while [[ "$1" > 0 ]]
do
  case $1 in
   -a| --assembly_file)
    shift
    ASSEMBLY=$1
    shift
  ;;
   -b| --alignment_file)
    shift
    ALIGNMENT=$1
    shift
  ;;
   -i| --indel_fraction)
    shift
    INDEL_FRAC=$1
    shift
  ;;
   -pt| --picard_path)
    shift
    PICARD_PATH=$1
    shift
  ;;
   -pl| --pilon_path)
    shift
    PILON_PATH=$1
    shift
  ;;
   -w| --warnings)
    shift
    WARNINGS=$1
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

if [[ -z "$WARNINGS" ]]; then

echo -e "\n"

echo "- The file $ASSEMBLY will be used as the reference assembly.

- The $ALIGNMENT will be used as alignment file.

- You selected $INDEL_FRAC as the indel fraction cut-off.

- pacbio_final_${PREFIX}.fasta will be the final corrected fasta."


##STARTING WITH PICARD...
echo -e "\n"
echo "Let's start with PICARD!!"

mkdir picard &>/dev/null
cd picard
ln -s ../$ASSEMBLY &>/dev/null
ln -s ../$ALIGNMENT &>/dev/null

#MarkDuplicates, AddOrReplaceReadGroups, and realigned processes.

java -Xmx10g -jar $PICARD_PATH MarkDuplicates INPUT=$ALIGNMENT OUTPUT=duplicate_${PREFIX}_PacBiovsIllumina.bam METRICS_FILE=duplicate_${PREFIX}_PacBiovsIllumina.bam.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 1>/dev/null 2>picard_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PICARD: MarkDuplicates!!For more information see picard_log.txt"
		exit 1
	fi 

java -Xmx10g -jar $PICARD_PATH AddOrReplaceReadGroups I=duplicate_${PREFIX}_PacBiovsIllumina.bam O=duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam RGID=[1,2] RGSM=1 RGLB=lib1 RGPU=unit1 RGPL=illumina 1>/dev/null 2>picard_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PICARD: AddOrReplaceReadGroups!!For more information see picard_log.txt"
		exit 1
	fi 

java -Xmx10g -jar $PICARD_PATH BuildBamIndex I=duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam 1>/dev/null 2>picard_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PICARD: BuildBamIndex!!For more information see picard_log.txt"
		exit 1
	fi 

rm -f duplicate_${PREFIX}_PacBiovsIllumina.bam duplicate_${PREFIX}_PacBiovsIllumina.bam.txt 1>/dev/null 2>picard_log.txt
find -type l -delete
cd ../

echo "
                                      _.._
                                     /,-. |_
                                 ___|/   '''
                           _____/___|______
                          /      ____      |
                         |      ('')))      |
                          |____() - -))____/
                               )( _>_)(
                              ()|'--()))
                             /'~~~/|~~'|
                            /   _/| )   |
                           /      /|     |
                         _/      |  |     |
                       _/        /  |      |
                      ///)|     (    )      )
                     / )(  )    )><><|      (
                     |_|(_,/    / ,  '|      )
                      ) | /    ( / | | )    /
                     (  | (    )_._._._(    |
                      ) |  |  /  /  |   |   |
                     (  |   |(   |  |   ||  |
                      ) |   _|__/    |__|_) |
                      |_|   |   /    |   /( |
                        |'. |_ /      | _| )|
                        |  ||- |      | -| |(
                        | ,--. |      | ,--.|
                    gnv | |____|''==--/____|-
"

echo "PICARD finished correctly!!"
#STARTING WITH PILON...

echo -e "\n"
echo "Continue with PILON!!"

mkdir pilon &>/dev/null
cd pilon
ln -s ../$ASSEMBLY &>/dev/null
ln -s ../picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam &>/dev/null

samtools index duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam
java -Xmx30g -jar $PILON_PATH --genome $ASSEMBLY --frags duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam  --output ${PREFIX}_pilon --changes --fix "indels" 1>/dev/null 2>pilon_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PILON!! For more information see pilon_log.txt"
		exit 1
	fi 

python3 ${ARAMIS_PATH}parser_pilon_bed.py ${PREFIX}_pilon.changes #parser pilon output to bed


	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE THE PARSER SCRIPT FROM PILON TO BED FILE!!"
		exit 1
	fi 


rm -f ${PREFIX}_pilon.fasta ${PREFIX}_pilon.changes  &>/dev/null
find -type l -delete
cd ..

echo "
                                                _..._
                                          ,--. /.---.|
                                     .---/____||--.  '
                                    (    '----'    )
                                     '-..........-'
                                       __|'--(__
                                      /'~~~/|~~'|
                                     /   _/| )   |
                                    /      /|     |
                                  _/      |  |     |
                                _/        /  |      |
                              _///)|     (    )      )
                             / )/   )    )><><|      (
                             L/(_,  /    / ,  '|      )
                              /    /    ( / | | )    /
                             /(    (    )_._._._(    |
                            /  )    |  /  /  |   |   |
                           /  (      |(   |  |   ||  |
                          /    )     _|__/    |__|_) |
                         /     |__   |   /    |   /( |
                        /         '. |_ /      | _| )|
                       /            ||- |      | -| |(
                      /            ,--. |      | ,--.||
                     /          gnv|____|''==--/____|-'
"


echo -e "\n"
echo "PILON finished correctly !!"

#STARTING WITH PACBIO-UTILITIES...

echo "Finally PACBIO-UTILITIES!!"

mkdir pacbio_utilities &>/dev/null
cd pacbio_utilities
ln -s ../$ASSEMBLY &>/dev/null
ln -s ../picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam &>/dev/null
ln -s ../pilon/${PREFIX}_pilon.changes.bed &>/dev/null

#Detecting all indels with specified indel fraction.
pacbio-util indel-targets -f $ASSEMBLY duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam --include-bad-indels --indel-frac $INDEL_FRAC > targets_All_${INDEL_FRAC}.txt 2>pacbioUtilities_log.txt #Bad indels included

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-targets!! For more information see pacbioUtilities_log.txt"
		exit 1
	fi 

#Generating files with common and not common indels position between pilon and pacbio_utilities.
awk -v infr=$INDEL_FRAC 'BEGIN{FS=OFS="\t"} $7>=infr' targets_All_${INDEL_FRAC}.txt | grep 'bad' > targets_onlybad_${INDEL_FRAC}.txt
awk -v infr=$INDEL_FRAC 'BEGIN{FS=OFS="\t"} $7>=infr' targets_All_${INDEL_FRAC}.txt | grep 'good' > targets_onlygood_${INDEL_FRAC}.txt

python3 ${ARAMIS_PATH}PilonCheck.py ${PREFIX}_pilon.changes.bed targets_onlybad_${INDEL_FRAC}.txt pilon_common.txt pilon_not_common.txt not_warning.txt warning.txt 2>pacbioUtilities_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE THE SCRIPT CREATED TO CHECK BAD PACBIO-UTILITIES IN PILON OUTPUT !!"
		exit 1
	fi 


sed '2,$d' targets_All_${INDEL_FRAC}.txt > header
cat targets_onlygood_${INDEL_FRAC}.txt pilon_common.txt | sort -k2,2n > targets_final_${PREFIX}.txt

#create list file with contigs names.
grep '>' $ASSEMBLY |sed 's/>//g' > list.txt

#Sort targets file following the fasta file order
FILE_TO_SORT="targets_final_${PREFIX}.txt"
LIST_FILE="list.txt"
TMP_FILE=$(mktemp)

while read LINE; do
    grep "$LINE" "$FILE_TO_SORT" |awk 'BEGIN{FS=OFS="\t"} $9 = toupper($9)' >> "$TMP_FILE"
done <"$LIST_FILE"

mv -f "$TMP_FILE" "$FILE_TO_SORT"

cat header targets_final_${PREFIX}.txt > targets_final_${PREFIX}_sorted.txt



echo "
                                       _..._
                                 ,::. /.---.|
                            .::-/___:||::.  '
                           (::  '----' .::)
                            '=:.......::='
                              __|'--(__
                             /::.~/|~::|
                            /::._/| ).::|
                           /::.   /|  .::|
                         _/::.   |  |  .::|
                       _/:::.    /  |   .::|
                      ///)|:.   (    )   .::)
                     / )(  ):.  )><><|   .::(
                     |_|(_,/:.  / ,  '|   .::)
                      ):|:/::. ( / | | )  .:/
                     (::|:(:::.)_._._._(  .:|
                      ):|::|::/  /::|   |..:|
                     (::|:::|(   |::|   ||::|
                      ):|::::|__/::::|__|:):|
                      |:|:::|   /::::|   /(:|
                        |'::|_ /::::::| _|:)|
                        |  ||- |::::::| -|:|(
                        | ,--. |::::::| ,--.||
                    gnv | |____|''==::/____|-'
"

#Generating correcting fasta.
/usr/local/bin/pacbio-util indel-apply -f $ASSEMBLY -t targets_final_${PREFIX}_sorted.txt --include-bad-indels > pacbio_final_${PREFIX}.fasta 2>pacbioUtilities_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-apply!! For more information see pacbioUtilities_log.txt"
		exit 1
	fi 

cp pacbio_final_${PREFIX}.fasta ../ &>/dev/null
cp targets_final_${PREFIX}_sorted.txt ../ &>/dev/null
find -type l -delete
rm -f header list.txt targets_final_${PREFIX}.txt targets_All_${INDEL_FRAC}.txt &>/dev/null
cd ../


echo -e "\n"
echo "PACBIO-UTILITILES finished correctly !!"


else

   echo "- The file $ASSEMBLY will be used as the reference assembly.
   
  - The file $WARNINGS will be used as the new file of targets to be corrected.

  - pacbio_final_withwarnings_${PREFIX}.fasta will be the final corrected fasta."
 
  ##STARTING WITH ADDITIONAL CORRECTION...
  echo -e "\n"
  echo "Let's do the new correction !!"

  mkdir second_correction &>/dev/null
  cd second_correction
  ln -s ../$ASSEMBLY &>/dev/null
  ln -s ../$WARNINGS &>/dev/null
  ln -s ../pacbio_utilities/targets_onlygood_${INDEL_FRAC}.txt &>/dev/null
  ln -s ../pacbio_utilities/not_warning.txt &>/dev/null
  
  echo "#assembly:$ASSEMBLY" > header
  cat header targets_onlygood_${INDEL_FRAC}.txt not_warning.txt $WARNINGS | sort -k2,2n > targets_final_withwarnings_${PREFIX}.txt
 
#Let's re-run pacbio-util  indel-apply

  pacbio-util indel-apply -f $ASSEMBLY -t targets_final_withwarnings_${PREFIX}.txt > pacbio_final_withwarnings_${PREFIX}.fasta

    if [[ $? != 0 ]]; then
        echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-apply!!"
        exit 1
    fi

  cp pacbio_final_withwarnings_${PREFIX}.fasta ../ &>/dev/null
  find -type l -delete

fi


echo "
                                  /   |
 _                        )      ((   ))     (
(@)                      /||      ))_((     /||
|-|                     / | |    (/||/|)   / | |                      (@)
| | -------------------/--|-voV---|'|'/--Vov-|--|---------------------|-|
|-|                         '^'   (o o)  '^'                          | |
| |                               '|Y/'                               |-|
|-|                                                                   | |
| |                         That's all folks!!!                       |-|
|-|                                                                   | |
| |  Eva Sacristán Horcajada (CBMSO)                                  |-|
| |            Sandra González de la Fuente (CBMSO)                   | |
|-|                         Ramón Peiró Pastor  (CBMSO)               |-|
|_|___________________________________________________________________| |
(@)              l   /| /         ( (       | /|   l                '||-|
                 l /   V           | |       V   | l                  (@)
                 l/                _) )_          |I
                                   '| /'
"
