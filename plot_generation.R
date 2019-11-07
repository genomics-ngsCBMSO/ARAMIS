#!/usr/bin/env Rscript

##CLEAN AND LIBRARIES


print(" Checking if all libraries are installed ")
rm(list=ls())

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

if ("optparse" %in% row.names(installed.packages())  == FALSE) install.packages("optparse")
if ("ggplot2" %in% row.names(installed.packages())  == FALSE) install.packages("ggplot2")
if ("cowplot" %in% row.names(installed.packages())  == FALSE) install.packages("cowplot")
if ("grid" %in% row.names(installed.packages())  == FALSE) install.packages("grid")
if ("epicontacts" %in% row.names(installed.packages())  == FALSE) install.packages("epicontacts")
if ("gridExtra" %in% row.names(installed.packages())  == FALSE) install.packages("gridExtra")
library("optparse")
library("ggplot2")
library("cowplot")
library("grid")
library("epicontacts")
library("gridExtra")

#SUBFUNCTIONS

snpposi.plot.indels <- function(x, genome.size, smooth=0.1, col="royalblue", alpha=.2)
{
        out <- ggplot(data.frame(x=x), aes(x=x)) + xlim(0, genome.size)
        out <- out + geom_density(adjust=smooth, fill=transp(col,alpha=alpha), colour=col) + geom_rug(colour=col,alpha=.7)
        out <- out + labs(x="Nucleotide position", title="Distribution of Indels")
    return(out)
}

#PARSE ARGUMENTS
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-l","--list"), type="character", default =NULL,
              help="table of chromosomes to analyse with their size", metavar="character"),
  make_option(c("-c", "--cov"), type="character", default="NULL", 
              help="file with coverage information separated by platform", metavar="character"),
  make_option(c("-g", "--gcskew"), type="character", default="NULL",
             help="gcskew of each chromosome",  metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("AT LEAST ONE ARGUMENT MUST BE PROVIDED", call.=FALSE)
}

## READING DATA

print(" Reading data ")

#ALL INFO FILE
info <- read.table(opt$file, sep="\t", header=T)

emptycols <- sapply(info, function (k) all(is.na(k)))

if (length(colnames(info[emptycols])) != 0) {
   print(paste("WARNING: ", colnames(info[emptycols]),"IS EMPTY!"))
   }

homtype <- as.data.frame(table(info$Type,info$Length))
colnames(homtype) <- c("Type","Length","Number")

size <- read.table(opt$list, sep="\t")

#COVERAGE FILE

cov <- read.table(opt$cov, header=T, sep="\t")

#GCSKEW FILE

gcs <- read.table(opt$gcskew, header=F, sep="\t")
colnames(gcs) <- c("Chromosome","Position","GC_Skew","Cumulative_GC_Skew")

## Generar los plots

print(" Generating plots ")

if (nrow(size) != 1) {
  G <- ggplot(homtype, aes(x=Length, y=Number, fill=Type)) + geom_bar(stat="identity", position=position_dodge()) 
  G <- G + labs(x="Length of the homopolymer", y = "Number of Homopolymers", coord_cartesian(ylim=c(0,50))) + theme_bw()
  ggsave("Variants_distribution_in_homopolymers_all.png", plot=G, device="png", width = 15, height = 8, dpi=352)
  ggsave("Variants_distribution_in_homopolymers_all.png", plot=G, device="pdf", width = 15, height = 8, dpi=352)
}


for (i in 1:nrow(size)) {
  genome.size <- size[i,2]
  var <- size[i,1]
  print(paste("Plots will be generated for", var))
  chro <-  subset(info, Chromosome == as.character(var))
  if (nrow(chro) == 0)
  {next}
  gcs_chro <- subset(gcs, Chromosome == as.character(var))
  if (nrow(gcs_chro) == 0)
  {next}
  #BARPLOT
  homtype_chro <- as.data.frame(table(chro$Type,chro$Length))
  if (nrow(homtype_chro) == 0)
  {next}
  colnames(homtype_chro) <- c("Type","Length","Number")
  G <- ggplot(homtype_chro, aes(x=Length, y=Number, fill=Type)) + geom_bar(stat="identity", position=position_dodge()) + labs(x="Length of the homopolymer", y = "Number of Homopolymers", title= paste("Distribution of variants in Homopolymers \n",var), coord_cartesian(ylim=c(0,50))) + theme_bw()
  ggsave(filename=paste("Variants_distribution_in_homopolymers_",var,".png", sep=""), plot=G, width = 15, height = 8, dpi=352)
  ggsave(filename=paste("Variants_distribution_in_homopolymers_",var,".pdf", sep=""), plot=G, width = 15, height = 8, dpi=352)

  #DENSITY PLOTS
  a <- snpposi.plot.indels(chro$Position, genome.size, smooth=0.1, col="royalblue", alpha=0.2)
  ggsave(filename=paste("Distribution_indels_",var,".png",sep=""), plot=a, device="png", width = 15, height = 8, dpi=352)
  ggsave(filename=paste("Distribution_indels_",var,".pdf", sep=""), plot=a, device="pdf", width = 15, height = 8, dpi=352)
  
  g_top <- ggplot(gcs_chro, aes(x=Position, y=Cumulative_GC_Skew)) + xlim(0, genome.size)+ geom_line() + labs(x="Nucleotide Position", y = "Cumulative_GC_Skew", title= paste("GC Skew of",var)) + theme_bw()
  together <- plot_grid(g_top, a, ncol = 1, nrow =2, labels=c("A","B"), axis = 'lb', rel_heights = c(1,2)) 
  ggsave(filename=paste("Distribution_indels_with_gcskew",var,".png",sep=""), plot=together, device="png", width = 15, height = 8, dpi=352)
  ggsave(filename=paste("Distribution_indels_with_gcskew",var,".pdf",sep=""), plot=together, device="pdf", width = 15, height = 8, dpi=352)

  ##DOT PLOT
  chro_cov <-  subset(cov, Chromosome == as.character(var))
  if (nrow(chro_cov) == 0)
  {next}
  out2 <- ggplot(chro_cov, aes(Position, Indel_fraction)) + xlim(0, genome.size) + geom_point(aes(colour=chro_cov$Platform)) 
  out2 <- out2 + labs(x="Nucleotide Position", y ="Indel Fraction", title =paste("Indel fraction by Sequencing Technology \n across", var), colour="Sequencing Technology") 
  out2 <- out2 + geom_hline(yintercept=as.numeric(0.5),col="darkgrey") + theme_bw()
  ggsave(filename=paste("Indel_fraction_",var,".png", sep=""), plot=out2, device="png", width = 15, height = 8, dpi=352)
  ggsave(filename=paste("Indel_fraction_",var,".pdf", sep=""), plot=out2, device="pdf", width = 15, height = 8, dpi=352)
}


##INTENTO PLOT
#nona <- subset(chro, !is.na(Type))
#out <- ggplot(nona, aes(x=Position, colour=Homopolymer)) + xlim(0, genome.size)
#out <- out + geom_density(adjust=0.1, size=0.8)
#out <- out + labs(x="Nucleotide position", title="") 
#ggsave(filename=paste("Indels_and_G4_lines",var,".png",sep=""), plot=out, device="png", width = 15, height = 8, dpi=352)


print(" All plots done!!! ")

