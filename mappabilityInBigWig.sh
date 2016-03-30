## Making the annotations files needed for mappability correction by the R package QDNAseq
## Jocelyn Sietsma Penington, Aug2015
## This is using the version of the reference genome with untransformed chromosome names 
cd /wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/malaria/vaccine/genomes/

gem-indexer -i PlasmoDB-12.0_Pfalciparum3D7_Genome.fasta -o Pfalciparum3D7_Genome

gem-mappability -I Pfalciparum3D7_Genome.gem -l 50 -o Pfalciparum3D7_Genome.l50default

gem-2-wig -I Pfalciparum3D7_Genome.gem -i Pfalciparum3D7_Genome.l50default.mappability \
 -o Pfalciparum3D7_Genome.l50default

## remove the portions of the chromosome headers of the form:
## | organism=Plasmodium_falciparum_3D7 | version=2013-03-01 | length=1200490 | SO=chromosome

sed 's/|.*|.*|.*|.*chromosome//' Pfalciparum3D7_Genome.l50default.wig > Pfalciparum3D7_Genome.l50default.wig2
sed 's/|.*|.*|.*|.*chromosome//' Pfalciparum3D7_Genome.l50default.sizes > Pfalciparum3D7_Genome.l50default.sizes2
mv Pfalciparum3D7_Genome.l50default.wig2 Pfalciparum3D7_Genome.l50default.wig
mv Pfalciparum3D7_Genome.l50default.sizes2 Pfalciparum3D7_Genome.l50default.sizes

wigToBigWig Pfalciparum3D7_Genome.l50default.wig Pfalciparum3D7_Genome.l50default.sizes Pfalciparum3D7_Genome50mer.bigWig

## Remove intermediate files
rm Pfalciparum3D7_Genome.gem
rm Pfalciparum3D7_Genome.l50default.*
