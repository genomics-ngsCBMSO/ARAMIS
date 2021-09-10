# ARAMIS

# Accurate long-Reads Assembly correction Method for Indel errorS

## INTRODUCTION

PacBio long reads can be used to generate de novo genome assemblies. However, we have noticed small insertions and deletions in those Pacbio-based assemblies, mostly single-base errors and frequently in homopolymer regions. 

Illumina reads on the other hand show a lower error rate and can be used to correct the possible single-based errors caused by PacBio long reads.

Thus, we have developed a collection of scripts to detect and correct those errors in PacBio assemblies using Illumina short reads. 
Here we described the pipeline we recommend in order to obtain a new corrected assembly in fasta format.  

 * For a full description of the module, visit the project page:
   Unpublished yet.

 * To submit bug reports and feature suggestions, or to track changes:
   genomicangs@cbm.csic.es

## CITATION

Thank you for using ARAMIS. Please, cite:

E Sacristán-Horcajada, S González-de la Fuente, R Peiró-Pastor, F Carrasco-Ramiro, R Amils, J M Requena, J Berenguer, B Aguado, ARAMIS: From systematic errors of NGS long reads to accurate assemblies, Briefings in Bioinformatics, 2021;, bbab170, https://doi.org/10.1093/bib/bbab170
  
## REQUIREMENTS

This pipeline requires the following dependencies:

* [Picard](https://github.com/broadinstitute/picard)

* [Pilon](https://github.com/broadinstitute/pilon)

* [samtools](https://github.com/samtools/samtools.git)

* [python3](https://www.python.org/download/releases/3.0/)

* [PacBio-utilities](https://github.com/douglasgscofield/PacBio-utilities)

* [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline)

* [seqkit](https://github.com/shenwei356/seqkit)

* gc_skew.py from [iRep](https://github.com/christophertbrown/iRep)

* fasta.py from [ctbBio](https://github.com/christophertbrown/bioscripts)

* [R](https://www.r-project.org/)

Picard, Pilon, PacBio-utilites, gc_skew.py and fasta.py are included in the folder 'dependencies' and do not need to be installed. Take into account that gc_skew.py was modified to use fasta.py from the dependencies folder.

The other dependencies can be installed inside a [conda](https://www.anaconda.com/) environment.

```
conda create -n Aramis python=3.8.5 r-base=4.0.5 scipy matplotlib seqkit samtools -c bioconda
conda activate Aramis
```

We recommend downloading [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline), decompress and add it to the path (include into .bashrc) because the conda package was in a lower version.

### R Packages

We recommend installing R packages manually by entering R through your terminal and type:

```
install.packages("optparse")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("grid")
install.packages("epicontacts")
install.packages("gridExtra")
```

## INSTALLATION

The downloaded folder contains 2 main scripts and auxiliary scripts needed to do the complete pipeline:

*Main scripts*

* correction.sh (AKA D'Artagnan)
* indel_analysis.sh (AKA Aramis)

*Auxiliary scripts*

* parser_pilon_bed.py (AKA Athos)
* PilonCheck.py (AKA Porthos)
* combine_info.py (AKA Richelieu)
* plot_generation.R (AKA Julieta)

ARAMIS's folder should be included in the path (we strongly recommend to add it permanently to /home/$USER/.bashrc) and the script must have execution permission:

```
export PATH="$PATH:/path/to/ARAMIS/" # Include in .bashrc
cd ARAMIS/
chmod -R +x ./
```

The script correction.sh assumes that the location of the .jar files of Picard Tools and Pilon, and PacBio-utilities is the folder dependencies. Please, do not erase this folder.

In case you want to use a different Picard or Pilon executable, add them as parameters when running the pipeline (--picard_path and --pilon_path)


## USAGE

The input files needed to run the pipeline are:

* **PacBio assembly fasta file:** De novo assembly obtained only with PacBio long reads in fasta format
* **Illumina alignment bam file:** Alignment of short Illumina reads against the PacBio assembly in bam format. We recommend using BWA aligner.
* **PacBio alignment bam file:** Alignment of long PacBio reads against the PacBio assembly in bam format. We recommend using Blasr aligner.

All of the input files for all the scripts MUST be in the work directory.  

1. Running the tests to correct the pacbio assembly file

```
correction.sh  -a PacBio.fasta -b Illumina.bam -i 0.5 -p prefix

OPTIONS
    -a  | --assembly_file   PacBio assembly in Fasta format (REQUIRED).
    -b  | --alignment_file  Input alignment in BAM (sorted) format (REQUIRED).
    -i  | --indel_fraction  Do not report indels at a position if the fraction of reads containing them is below FLOAT.
    -pt | --picard_path     To indicate a different path for picard (default: dependencies folder).
    -pl | --pilon_path      To indicate a different path for picard (default: dependencies folder).
    -w  | --warnings        Alternative pipeline after manual correction of indels flags as warning.
    -p  | --prefix          Prefix for output files (default: standard).
    -f  | --force           Option to force recomputation from the beginning (default: Inactive). Must be added at the end.
    -h  | --help            Display help.
    -c  | --citation        Display citation.
```

This step will generate 3 folders:

* picard/
    - duplicate_prefix_PacBiovsIllumina_readGroup.bam: Output of the pre-proccesing of the alignment
    - duplicate_prefix_PacBiovsIllumina_readGroup.bai: Index of bam file. Neccesary if you want to visualize the files in Genome Browers as IGV

* pilon/
    - prefix.pilon.changes.bed: Tabular file including the indel errors in the PacBio Assembly detected by PILON software

* pacbio_utilities/
    - pacbio_final_prefix.fasta: Final corrected genome aasembly in fasta format 
    - targets_final_sorted_prefix.txt: Final targets that passed the indel fraction threshold and were corrected in the genome assembly
    - warnings.txt: List of positions where indels were detected and corrected, but marked as warnings for possible second correction after manual checking
    - not_warning.txt: List of positions where indels were detected and corrected. Not necessary manual correction
    - targets_only_good.txt: List of positions flagged as good by pacbio-utilities directly
    - targets_only_bad.txt: List of positions flagged as bad by pacbio-utilities and thus check with Pilon



2. Running a second additional correction after modifiying the generate warnings file (OPTIONAL)


```
correction.sh -a PacBio.fasta -b Illumina.bam -i 0.5 -w warnings.txt -p prefix
```

This step will generate 1 folder:

   - pacbio_final_withwarnings_prefix.fasta: New genome assembly in fasta format generated after the second correction   



3. Running the tests to generate coverage, gc-skew and homopolymers statistics with their plots

```
indel_analysis.sh -a PacBio.fasta -b1 Illumina.bam -b2 PacBio.bam -p prefix -t <targets file>

OPTIONS
     -a | --assembly_file             PacBio assembly in Fasta format (REQUIRED).
    -b1 | --illumina_alignment_file   input illumina alignment in BAM (sorted) format (REQUIRED).
    -b2 | --pacbio_alignment_file     input pacbio alignment in BAM (sorted) format (REQUIRED).
     -p | --prefix                    Prefix for output files (default: standard)
     -t | --targets_file              Target file generated by correction.sh (REQUIRED).
     -h | --help                      Display help.
     -c | --citation                  Display citation.    
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


## MAINTAINERS

Current maintainers:

 * Eva Sacristán Horcajada (CBMSO) - esacristan@cbm.csic.es
 * Sandra González de la Fuente (CBMSO) - sandra.g@cbm.csic.es
 * Adrián Gómez Repollés (CBMSO) - adrian.gomez@cbm.csic.es

 ## Authors

 * Eva Sacristán Horcajada (CBMSO) - esacristan@cbm.csic.es
 * Sandra González de la Fuente (CBMSO) - sandra.g@cbm.csic.es
 * Ramón Peiró Pastor (CBMSO) - rpeiro@cbm.csic.es

This pipeline was developed by members of the Genomics and NGS Core Facility of the Centro de Biología Molecular Severo-Ochoa (CBMSO, Madrid, Spain),
a joint Center between the Consejo Superior de Investigaciones Científicas (CSIC) (Spanish Research Council) and the Universidad Autónoma de Madrid (UAM). 


## Acknowledgments

The authors appreciate the advice and helpful comments received from members of the Genomics and NGS Core Facility (GENGS, http://genomics-ngs.cbm.uam.es/) at the Centro de Biología Molecular Severo Ochoa (CBMSO, CSIC-UAM) which is part of the CEI UAM+CSIC, Madrid, Spain. Genome sequence data from European Nucleotide Archive (ENA) were invaluable for this work and their provision in the public domain is gratefully acknowledged. The authors thank all the staff members of Dr Berenguer, Dr Requena and Dr Amils laboratories at the CBMSO, especially Dr Alba Blesa, Dr Mercedes Sanchez, Esther Camacho and Jose Manuel Martínez. Computational time from the Centro de Computación Científica (CCC) of Universidad Autónoma de Madrid is also gratefully acknowledged.


