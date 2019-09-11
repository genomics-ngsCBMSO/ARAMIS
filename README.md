# ARAMIS

# Accurate long Reads Assemblies Method to correct Indel errors with Short reads.


## INTRODUCTION

PacBio long reads can be used to generate de novo genome assemblies. However, we have noticed small insertions
and deletions in those Pacbio-based assemblies, mostly single-base errors and frequently in homopolymer regions. 

Illumina reads on the other hand show a lower error rate and can be used to correct the possible single-based errors caused by PacBio long reads.

Thus, we have developed a collection of scripts to detect and correct those errors in PacBio assemblies using Illumina short reads. 
Here we described the pipeline we recommend in order to obtain a new corrected assembly in fasta format.  

 * For a full description of the module, visit the project page:
   github
   paper

 * To submit bug reports and feature suggestions, or to track changes:
   genomicangs@cbm.csic.es

## CITATION


  
## REQUIREMENTS

This pipeline requires the following dependences:

* [Picard](https://github.com/broadinstitute/picard)

* [Pilon](https://github.com/broadinstitute/pilon)

* [samtools](https://github.com/samtools/samtools.git)

* [python3](https://www.python.org/download/releases/3.0/)

* [PacBio-utilities](https://github.com/douglasgscofield/PacBio-utilities)

* [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline)

* [seqkit](https://github.com/shenwei356/seqkit)

* gc_skew.py from [iRep](https://github.com/christophertbrown/iRep)

* [R](https://www.r-project.org/)

### R Packages

We recommend installing R packages manually.

Enter R through your terminal and type:

```
install.packages("optparse")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("grid")
install.packages("epicontacts")
install.packages("gridExtra")

```

## INSTALLATION

Make sure you have installed all the dependencies before try to run the pipeline!

The downloaded folder contains 2 main scripts and a folder with all the additional scripts needed to do the complete pipeline:

*Main programs*

* correction.sh (AKA D'Artagnan)
* indel_analysis.sh (AKA Aramis)

*scripts/ folder*

* parser_pilon_bed.py (AKA Athos)
* PilonCheck.py (AKA Porthos)
* combine_info.py (AKA Richelieu)
* plot_generation.R (AKA Julieta)

All of the scripts must be located in /usr/local/bin and with execution permission:

```
cd ARAMIS/
chmod -R +x ./
sudo cp correction.sh indel_analysis.sh scripts/parser_pilon_bed.py scripts/PilonCheck.py scripts/combine_info.py scripts/plot_generation.R /usr/local/bin

```


The script correction.sh assumes that the path to the executables .jar files of Picard Tools and Pilon is:

/home/$USER/picard.jar (or ~/picard.jar)
/home/$USER/pilon.jar (or ~/pilon.jar)

If this is not true, modify the correspondant lines in the script indicating the correct path for each software, 
or add them as parameters when running the pipeline (--picard_path and --pilon_path)


## USAGE

The input files needed to run the pipeline are:

* *PacBio assembly fasta file:* De novo assembly obtained only with PacBio long reads in fasta format
* *Illumina alignment bam file:* Alignment of short Illumina reads against the PacBio assembly in bam format. We recommend using BWA aligner.
* *PacBio alignment bam file:* Alignment of long PacBio reads against the PacBio assembly in bam format. We recommend using Blasr aligner.

All of the input files for all the scripts MUST be in the work directory.  

1. Running the tests to correct the pacbio assembly file

```
correction.sh  -a PacBio.fasta -b Illumina.bam -i 0.5 -p prefix

OPTIONS
    -a  | --assembly_file     PacBio assembly in Fasta format (REQUIRED)
    -b  | --alignment_file    Input alignment in BAM (sorted) format (REQUIRED)
    -i  | --indel_fraction    Do not report indels at a position if the fraction of reads containing them is below this number FLOAT
    -pt | --picard_path       If picard.jar is not in the path /home/$USER , please indicate the correct path with this option 
    -pl | --pilon_path        If pilon.jar is not in the path /home/$USER , please indicate the correct path with this option 
    -w  | --warnings          Alternative pipeline after manual correction of indels flags as warning
    -p  | --prefix            Prefix for output files
    -h  | --help              Display help

```

This step will generate 3 folders:

* picard
    - duplicate_prefix_PacBiovsIllumina_readGroup.bam: Output of the pre-proccesing of the alignment
    - duplicate_prefix_PacBiovsIllumina_readGroup.bai: Index of bam file. Neccesary if you want to visualize the files in Genome Browers as IGV

* pilon
    - prefix.pilon.changes.bed: Tabular file including the indel errors in the PacBio Assembly detected by PILON software

* pacbio_utilities
    - pacbio_final_prefix.fasta: Final corrected genome aasembly in fasta format 
    - targets_final_sorted_prefix.txt: Final targets that passed the indel fraction threshold and were corrected in the genome assembly
    - manual_correction_pilon_not_common_indelfraction.txt: List of positions where indels were detected but not corrected as they didn't pass all the filters
    - warnings.txt: List of positions where indels were detected and corrected, but marked as warnings for possible second correction after manual checking
    - not_warning.txt: List of positions where indels were detected and corrected. Not necessary manual correction
    - targets_only_good.txt: List of positions flagged as good by pacbio-utilities directly
    - targets_only_bad.txt: List of positions flagged as bad by pacbio-utilities and thus check with Pilon



2. Running a second additional correction after modifiying the generate warnings file (OPTIONAL)

```
correction.sh -a PacBio.fasta -b Illumina.bam -i 0.5 -w warnings.txt -p prefix
```

This step will generate 1 folder

    - pacbio_final_withwarnings_prefix.fasta: New genome assembly in fasta format generated after the second correction   



3. Running the tests to generate coverage, gc-skew and homopolymers statistics with their plots

```
indel_analysis.sh -a <PacBio assembly fasta file> -b1 <Illumina alignment bam file> -b2 <PacBio alignment bam file> -p <prefix for output files> -t <targets file>

OPTIONS
    -a  | --assembly_file             PacBio assembly in Fasta format
    -b1 | --illumina_alignment_file   Input Illumina alignment in BAM (sorted) format
    -b2 | --pacbio_alignment_file     Input PacBio alignment in BAM (sorted) format
    -p  | --prefix                    Prefix for output files
    -t  | --targets_file              Target file generated by correction.sh
    -h  | --help                      Display help
    
```

This step will generate 1 folder

* indel_information
    - prefix_homopolymers.txt: Position and information of all the homopolymers of at least 2 Cs,Gs,As or Ts in the genome
    - Variants_distribution_in_homopolymers_chromosome.png: Bar plot representing the homopolymers affected by the indels ordered by length and base
    - prefix_Illumina_coverage.wig: Coverage calculation of the Illumina alignment
    - prefix_PacBio_coverage.wig: Coverage calculation of the PacBio alignment
    - prefix_gc_skew.txt: Gc_skew calculation of the non-corrected genome assembly
    - Distribution_indels_cromosoma.pdf/png: Density plot representing the distribution of indels across the genome assembly.
    - Distribution_indels_withgcskew_cromosoma.pdf/png: Same plot as the last one but with gc skew information add
    - stats.txt: Stats calculated from the homopolymers information.
    - Indel_fraction_cromosoma.pdf: Plot representing the indel fraction of the Illumina alignment and the PacBio alignment 
    - prefix_allinfo.txt: Tabular file including all the information calculated for each position where an indel was detected and corrected

```

## MAINTAINERS

Current maintainers:
 * Eva Sacristán Horcajada (CBMSO) - esacristan@cbm.csic.es
 * Sandra González de la Fuente (CBMSO) - sandra.g@cbm.csic.es

 ## Authors

 * Eva Sacristán Horcajada (CBMSO) - esacristan@cbm.csic.es
 * Sandra González de la Fuente (CBMSO) - sandra.g@cbm.csic.es
 * Ramón Peiró Pastor (CBMSO) - rpeiro@cbm.csic.es

This pipeline was developed by memebers of the Genomics and NGS Core Facility of the Centro de Biología Molecular Severo-Ochoa (CBMSO, Madrid, Spain),
a joint Center between the Consejo Superior de Investigaciones Científicas (CSIC) (Spanish Research Council) and the Universidad Autónoma de Madrid (UAM). 


## Acknowledgments

* Agradecimientos del paper
