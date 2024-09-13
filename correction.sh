#!/bin/bash

function usage()
{
    echo "Usage: $0 -a <ASSEMBLY_FILE> -b <ALIGNMENT_FILE> -i <INDEL_FRAC> -p <PREFIX>"
    printf "\n"
    echo "In case you want to perform manual correction of warning targets, please run Usage2."
    echo "Usage2: $0 -a <ASSEMBLY_FILE> -w  <WARNINGS> -p <PREFIX>"
    printf "\n"
    echo "  -a  | --assembly_file   PacBio assembly in Fasta format (REQUIRED)."
    echo "  -b  | --alignment_file  Input alignment in BAM (sorted) format (REQUIRED)."
    echo "  -i  | --indel_fraction  Do not report indels at a position if the fraction of reads containing them is below FLOAT."
    echo "  -pt | --picard_path     To indicate a different path for picard (default: dependencies folder)."  
    echo "  -pl | --pilon_path      To indicate a different path for picard (default: dependencies folder)." 
    echo "  -w  | --warnings        Alternative pipeline after manual correction of indels flags as warning."
    echo "  -p  | --prefix          Prefix for output files (default: standard)."
    echo "  -f  | --force           Option to force recomputation from the beginning (default: Inactive). Must be added at the end."
    echo "  -h  | --help            Display help."
    echo "  -c  | --citation        Display citation."
}

function citation()
{
  echo "Thank you for using ARAMIS. Please, cite:"
  echo ""
  echo "E Sacristán-Horcajada, S González-de la Fuente, R Peiró-Pastor, F Carrasco-Ramiro, R Amils, J M Requena, J Berenguer, B Aguado, ARAMIS: From systematic errors of NGS long reads to accurate assemblies, Briefings in Bioinformatics, 2021;, bbab170, https://doi.org/10.1093/bib/bbab170"
  echo ""
}

# Default variables
ASSEMBLY=''
ALIGNMENT=''
INDEL_FRAC=0.5
WARNINGS=''
PREFIX='standard'
FORCE=false

# display usage if 
if [[ $# -lt 4  && $1 != "-h" && $1 != "--help" && $1 != "-c" && $1 != "--citation" ]]; then
  echo "ERROR: You must provided at least the assembly and alignment files."
  usage
  exit 2
fi

#Get parameters
while [[ "$1" > 0 ]]; do
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
    -p| --prefix)
      shift
      PREFIX=$1
      shift
      ;;
    -f| --force)
      shift
      FORCE=true
      shift
      ;;
    -h| --help)
      usage
      exit
      ;;
    -c | --citation)
      citation
      exit
      ;;
    *)
      echo "ERROR: Missing parameter!!!"
      exit
      ;;
  esac
done

function clear_links()
{
# Get link name from full path
assembly_link=$(basename $ASSEMBLY)
alingment_link=$(basename $ALIGNMENT)
[ ! -z "$WARNINGS" ] && warnings_link=$(basename $WARNINGS)
# Clean links
links=($assembly_link $alingment_link duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam ${PREFIX}_pilon.changes.bed $warnings_link targets_onlygood_${INDEL_FRAC}.txt not_warning.txt)
  for lnk in ${links[@]}
  do
    if [ -h "$lnk" ]
    then
      unlink $lnk
    fi
  done
}

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
  printf "\nExecution halted by user.\n"
  clear_links
  exit
}

# Get ARAMIS path
ARAMIS_PATH=$(/usr/bin/dirname $(/usr/bin/realpath $0))

# Check software dependencies
software=(pacbio-util pilon-1.24.jar picard.jar)
for pckg in ${software[@]}
do 
  if [[ ! -f "${ARAMIS_PATH}/dependencies/${pckg}" ]]
  then
    printf "\nERROR: dependencies missing (${pckg}).\n"
    exit
  else
    if [[ ${pckg} == 'pacbio-util' ]]
    then
      PBUTIL_PATH="${ARAMIS_PATH}/dependencies/${pckg}"
    elif [[ ${pckg} == 'picard.jar' && -z $PICARD_PATH ]]
    then
      PICARD_PATH="${ARAMIS_PATH}/dependencies/${pckg}"
    elif [[ ${pckg} == 'pilon-1.24.jar' && -z $PILON_PATH ]]
    then
      PILON_PATH="${ARAMIS_PATH}/dependencies/${pckg}"
    fi
  fi
done

# Check Samtools
type samtools 2> /dev/null 1>&2 
if [ $? != 0 ]
then
  printf "\nERROR: samtools missing. Install the required dependencies. We recommend installing a conda environment.\n"
  exit
else
  SAMTOOLS=$(type -p samtools)
fi

# Check Perl
type perl 2> /dev/null 1>&2
if [ $? != 0 ]
then
  printf "\nERROR: Perl missing. Install the required dependencies."
  exit
fi

# Check BioPerl
perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION' &> /dev/null
if [ $? != 0 ]
then
  printf "\nERROR: BioPerl missing. Install the required dependencies. We recommend installing a conda environment.\n"
  exit
fi

# Check if computation was already done
run_picard=true
run_pilon=true
run_pbutil=true
run_pbutil2=true
[[ -d picard && ${FORCE} == 'false' ]] && run_picard=false
[[ -d pilon && ${FORCE} == 'false' ]] && run_pilon=false
[[ -d pacbio_utilities && ${FORCE} == 'false' ]] && run_pbutil=false
[[ -d second_correction && ${FORCE} == 'false' ]] && run_pbutil2=false

# Start computation
if [[ -z "$WARNINGS" ]]
then

  printf "\n"
  echo "  -The file $ASSEMBLY will be used as the reference assembly.
  -The $ALIGNMENT will be used as alignment file.
  -You selected $INDEL_FRAC as the indel fraction cut-off.
  -pacbio_final_${PREFIX}.fasta will be the final corrected fasta."

##STARTING WITH PICARD...
if [[ ${run_picard} != true ]]
then
  printf "\n"
  echo "Warning: Picard already computed. Use the option '--force' to recompute the full correction or remove picard folder to compute again."
  printf "\n"
else
  printf "\n"
  echo "Let's start with PICARD!!"

  mkdir picard &>/dev/null
  # cd picard
  ln -s $ASSEMBLY &>/dev/null
  ln -s $ALIGNMENT &>/dev/null

  #MarkDuplicates, AddOrReplaceReadGroups, and realigned processes.
  java -Xmx10g -jar $PICARD_PATH MarkDuplicates INPUT=$ALIGNMENT OUTPUT=picard/duplicate_${PREFIX}_PacBiovsIllumina.bam METRICS_FILE=picard/duplicate_${PREFIX}_PacBiovsIllumina.bam.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 1>/dev/null 2>picard/picard_log.txt

  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE PICARD: MarkDuplicates!!For more information see picard_log.txt"
  		exit 1
  	fi 

  java -Xmx10g -jar $PICARD_PATH AddOrReplaceReadGroups I=picard/duplicate_${PREFIX}_PacBiovsIllumina.bam O=picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam RGID=[1,2] RGSM=1 RGLB=lib1 RGPU=unit1 RGPL=illumina 1>/dev/null 2>picard/picard_log.txt

  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE PICARD: AddOrReplaceReadGroups!!For more information see picard_log.txt"
  		exit 1
  	fi 

  java -Xmx10g -jar $PICARD_PATH BuildBamIndex I=picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam 1>/dev/null 2>picard/picard_log.txt

  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE PICARD: BuildBamIndex!!For more information see picard_log.txt"
  		exit 1
  	fi 

  rm -f picard/duplicate_${PREFIX}_PacBiovsIllumina.bam picard/duplicate_${PREFIX}_PacBiovsIllumina.bam.txt 1>/dev/null 2>picard/picard_log.txt
  clear_links

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
  printf "\n"
  echo "PICARD finished correctly!!"
fi

#STARTING WITH PILON...
if [[ ${run_pilon} != true ]]
then
  echo "Warning: Pilon already computed. Use the option '--force' to recompute the full correction or remove pilon folder to compute again."
  printf "\n"
else
  printf "\n"
  echo "Continue with PILON!!"

  mkdir pilon &>/dev/null
  ln -s $ASSEMBLY &>/dev/null
  ln -s picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam &>/dev/null

  $SAMTOOLS index picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam
  java -Xmx30g -jar $PILON_PATH --genome $ASSEMBLY --frags picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam  --output pilon/${PREFIX}_pilon --changes --fix "indels" 1>/dev/null 2>pilon/pilon_log.txt

  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE PILON!! For more information see pilon_log.txt"
  		exit 1
  	fi 

  python3 ${ARAMIS_PATH}/parser_pilon_bed.py pilon/${PREFIX}_pilon.changes #parser pilon output to bed


  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE THE PARSER SCRIPT FROM PILON TO BED FILE!!"
  		exit 1
  	fi 

  rm -f pilon/${PREFIX}_pilon.fasta pilon/${PREFIX}_pilon.changes  &>/dev/null
  clear_links

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


  printf "\n"
  echo "PILON finished correctly !!"
fi

# STARTING WITH PACBIO-UTILITIES...
if [[ ${run_pbutil} != true ]]
then
  echo "Warning: PacBio-Utilities already computed. Use the option '--force' to recompute the full correction or remove pacbio_utilities folder to compute again."
  printf "\n"
  exit
else
  printf "\n"
  echo "Finally PACBIO-UTILITIES!!"

  mkdir pacbio_utilities &>/dev/null
  #cd pacbio_utilities
  ln -s $ASSEMBLY &>/dev/null
  ln -s picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam &>/dev/null
  ln -s pilon/${PREFIX}_pilon.changes.bed &>/dev/null

  # Detecting all indels with specified indel fraction.
  ${PBUTIL_PATH} indel-targets -f $ASSEMBLY picard/duplicate_${PREFIX}_PacBiovsIllumina_readGroup.bam --include-bad-indels --indel-frac $INDEL_FRAC > pacbio_utilities/targets_All_${INDEL_FRAC}.txt 2>pacbio_utilities/pacbioUtilities_log.txt #Bad indels included
  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-targets!! For more information see pacbioUtilities_log.txt"
  		exit 1
  	fi 

  #Generating files with common and not common indels position between pilon and pacbio_utilities.
  awk -v infr=$INDEL_FRAC 'BEGIN{FS=OFS="\t"} $7>=infr' pacbio_utilities/targets_All_${INDEL_FRAC}.txt | grep 'bad' > pacbio_utilities/targets_onlybad_${INDEL_FRAC}.txt
  awk -v infr=$INDEL_FRAC 'BEGIN{FS=OFS="\t"} $7>=infr' pacbio_utilities/targets_All_${INDEL_FRAC}.txt | grep 'good' > pacbio_utilities/targets_onlygood_${INDEL_FRAC}.txt

  python3 ${ARAMIS_PATH}/PilonCheck.py ${PREFIX}_pilon.changes.bed pacbio_utilities/targets_onlybad_${INDEL_FRAC}.txt pacbio_utilities/pilon_common.txt pacbio_utilities/pilon_not_common.txt pacbio_utilities/not_warning.txt pacbio_utilities/warning.txt 2>pacbio_utilities/pacbioUtilities_log.txt

  	if [[ $? != 0 ]]; then 
  		echo "ERROR TRYING TO EXECUTE THE SCRIPT CREATED TO CHECK BAD PACBIO-UTILITIES IN PILON OUTPUT !!"
  		exit 1
  	fi 


  sed '2,$d' pacbio_utilities/targets_All_${INDEL_FRAC}.txt > pacbio_utilities/header
  cat pacbio_utilities/targets_onlygood_${INDEL_FRAC}.txt pacbio_utilities/pilon_common.txt | sort -k1,1 -k2,2n > pacbio_utilities/targets_final_${PREFIX}.txt

  # Create list file with contigs names.
  grep '>' $ASSEMBLY |sed 's/>//g' > pacbio_utilities/list.txt

  # Sort targets file following the fasta file order
  FILE_TO_SORT="pacbio_utilities/targets_final_${PREFIX}.txt"
  LIST_FILE="pacbio_utilities/list.txt"
  TMP_FILE=$(mktemp)

  awk 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1], $0;}' $LIST_FILE $FILE_TO_SORT | sort -k 1,1n -k 3,3n | cut -d' ' -f2- | awk 'BEGIN{FS=OFS="\t"} $9 = toupper($9)' >> $TMP_FILE

  mv $TMP_FILE $FILE_TO_SORT

  cat pacbio_utilities/header pacbio_utilities/targets_final_${PREFIX}.txt > pacbio_utilities/targets_final_${PREFIX}_sorted.txt
fi

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

# Generate corrected fasta.
${PBUTIL_PATH} indel-apply -f $ASSEMBLY -t pacbio_utilities/targets_final_${PREFIX}_sorted.txt --include-bad-indels > pacbio_utilities/pacbio_final_${PREFIX}.fasta 2>pacbio_utilities/pacbioUtilities_log.txt

	if [[ $? != 0 ]]; then 
		echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-apply!! For more information see pacbioUtilities_log.txt"
		exit 1
	fi 

clear_links
rm -f pacbio_utilities/header pacbio_utilities/list.txt pacbio_utilities/targets_final_${PREFIX}.txt pacbio_utilities/targets_All_${INDEL_FRAC}.txt &>/dev/null

printf "\n"
echo "PACBIO-UTILITILES finished correctly!!"

else
  if [[ ${run_pbutil2} != true ]]
  then
    printf "\n"
    echo "Second correction already computed. Use the option '--force' to recompute the second correction or remove second_correction folder to compute again."
    exit
  else
    # Removing folder
    if [[ -d second_correction ]]
    then
      rm -rf second_correction
    fi

    printf "\n"
    echo "    -The file $ASSEMBLY will be used as the reference assembly. 
    -The file $WARNINGS will be used as the new file of targets to be corrected.
    -pacbio_final_withwarnings_${PREFIX}.fasta will be the final corrected fasta."

    ##STARTING WITH ADDITIONAL CORRECTION...
    printf "\n"
    echo "Let's do the new correction!!"

    mkdir -p second_correction &>/dev/null
    #cd second_correction
    ln -s $ASSEMBLY &>/dev/null
    ln -s $WARNINGS &>/dev/null
    ln -s pacbio_utilities/targets_onlygood_${INDEL_FRAC}.txt &>/dev/null
    ln -s pacbio_utilities/not_warning.txt &>/dev/null
    
    grep '>' $ASSEMBLY |sed 's/>//g' > second_correction/list.txt
    LIST_FILE="second_correction/list.txt"
    FILE_TO_SORT="second_correction/targets_final_withwarnings_${PREFIX}.txt"
    TMP_FILE=$(mktemp)

    echo "#assembly:$ASSEMBLY" > second_correction/header
    cat targets_onlygood_${INDEL_FRAC}.txt not_warning.txt $WARNINGS > $TMP_FILE
    awk 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1], $0;}' $LIST_FILE $TMP_FILE | sort -k 1,1n -k 3,3n | cut -d' ' -f2- | awk 'BEGIN{FS=OFS="\t"} $9 = toupper($9)' >> $FILE_TO_SORT
    cat <(cat second_correction/header) $FILE_TO_SORT > $TMP_FILE # Add header to sort file
    mv $TMP_FILE $FILE_TO_SORT

    #Let's re-run pacbio-util  indel-apply
    ${PBUTIL_PATH} indel-apply -f $ASSEMBLY -t second_correction/targets_final_withwarnings_${PREFIX}.txt > second_correction/pacbio_final_withwarnings_${PREFIX}.fasta

      if [[ $? != 0 ]]; then
          echo "ERROR TRYING TO EXECUTE PACBIO-UTILITIES: indel-apply!!"
          exit 1
      fi

    rm second_correction/list.txt
    clear_links
  fi
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
