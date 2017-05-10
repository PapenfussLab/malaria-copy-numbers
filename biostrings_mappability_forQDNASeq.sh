## Making the files needed to use QDNASeq R package with release 29 of Plasmodium
## genome from PlasmoDB

## Jocelyn Sietsma Penington 
## February 2017

PAPENDIR=/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/
REFDIR=$PAPENDIR/reference_genomes/plasmodium/PlasmoDB-29_Pfalciparum3D7

cd $REFDIR
## forgeBSgenomeDataPkg wants input in twoBit format
faToTwoBit PlasmoDB-29_Pfalciparum3D7_Genome.fasta   \
           PlasmoDB-29_Pfalciparum3D7_Genome.2bit
           
## Mappability
gem-indexer -i PlasmoDB-29_Pfalciparum3D7_Genome.fasta -o Pfalciparum3D7_Genome
gem-mappability -I Pfalciparum3D7_Genome.gem -l 50 -o Pfalciparum3D7_Genome.l50default
gem-2-wig -I Pfalciparum3D7_Genome.gem -i Pfalciparum3D7_Genome.l50default.mappability \
  -o Pfalciparum3D7_Genome.l50default

## remove the portions of the chromosome headers of the form:
## | organism=Plasmodium_falciparum_3D7 | version=2013-03-01 | length=1200490 | SO=chromosome
## because wigToBigWig fails on reading '|'

sed 's/|.*|.*|.*|.*chromosome//' Pfalciparum3D7_Genome.l50default.wig > Pfalciparum3D7_Genome.l50default.wig2
sed 's/|.*|.*|.*|.*chromosome//' Pfalciparum3D7_Genome.l50default.sizes > Pfalciparum3D7_Genome.l50default.sizes2
mv Pfalciparum3D7_Genome.l50default.wig2 Pfalciparum3D7_Genome.l50default.wig
mv Pfalciparum3D7_Genome.l50default.sizes2 Pfalciparum3D7_Genome.l50default.sizes

wigToBigWig Pfalciparum3D7_Genome.l50default.wig Pfalciparum3D7_Genome.l50default.sizes Pfalciparum3D7_Genome50mer.bigWig

## Remove intermediate files
if [[ -f Pfalciparum3D7_Genome50mer.bigWig ]] ; then
   rm Pfalciparum3D7_Genome.gem
   rm Pfalciparum3D7_Genome.l50default.*
fi

## Text file BioStringsFormatSeedFile.txt has been put in REFDIR
## forgeBSgenomeDataPkg("BioStringsFormatSeedFile.txt") is done in R  
## Following steps can be done in R using Wickham's devtools(), so no longer used here   

# R CMD build BSgenome.Pfalciparum3D7.PlasmoDB_29
# R CMD check BSgenome.Pfalciparum3D7.PlasmoDB_29 
# R CMD INSTALL -l /home/users/allstaff/penington.j/R/x86_64-pc-linux-gnu-library/3.3/  \
#    BSgenome.Pfalciparum3D7.PlasmoDB_29_0.1-0.tar.gz 
