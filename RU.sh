#!/bin/bash

mkdir ~/workdir/assignmet && cd ~/workdir/assignmet
work_dir="$(pwd)"

#------------------------------------------------------------
# 1- Download and Index the reference human genome
mkdir genome-data && cd genome-data
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz

# indexing
mkdir idx/
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir idx/ --genomeFastaFiles GCA_000001405.28_GRCh38.p13_genomic.fna

index="$work_dir/genome-data/idx/"

#------------------------------------------------------------
# 2- Download the samples and Prepare sample info.
#wget https://eu-gra.uploadfiles.io/get/kc0qqbvd
unzip ngs2-assignment-data.zip

R1="$work_dir/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq.gz"
R2="$work_dir/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq.gz"

SM="SRR8797509"
PL="Illumina"



#------------------------------------------------------------
# 3- Map the sample to the ref.



#------------------------------------------------------------
# 4- Mark duplicates and sort (Picard tool)


#------------------------------------------------------------
# 5- 
