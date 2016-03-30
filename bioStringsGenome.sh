faToTwoBit PlasmoDB-12.0_Pfalciparum3D7_Genome.fasta PlasmoDB-12.0_Pfalciparum3D7_Genome.2bit

vi BioStringsFormatSeedFile.txt 
# Make any necessary changes. Check filenames, chromosome names in circ_seqs field

## In R:
# library(BSgenome)
# setwd("/wehisan/home/allstaff/p/penington.j/Pfvaccine_ln/genomes") 
# forgeBSgenomeDataPkg("BioStringsFormatSeedFile.txt") 

R CMD build BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3
R CMD check BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3
R CMD INSTALL BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3_0.1-0.tar.gz 

## In R:
# library("BSgenome.Pfalciparum3D7.PlasmoDB.3D7v3")
