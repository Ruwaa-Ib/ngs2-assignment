#!/bin/bash

cd ~/workdir
git clone https://github.com/Ruwaa-Ib/ngs2-assignment.git
cd ngs2-assignment
work_dir="$(pwd)"

#------------------------------------------------------------
# 1- Download and Index the reference human genome
mkdir genome-data && cd genome-data

# download ref.
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz

# unzipping
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz
ref_fna="$work_dir/genome-data/GCA_000001405.28_GRCh38.p13_genomic.fna"

# hisat indexing
mkdir idx
hisat2-build $ref_fna idx/grch38
index="$work_dir/genome-data/idx/grch38"

# create a ref. index with faidx (.fai)
samtools faidx $ref_fna

# create a ref. dict. with Picard (.dict)
picard_path=$CONDA_PREFIX/share/picard-2.19.2-0
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=$ref_fna O=$ref_fna.dict

#------------------------------------------------------------
# 2- Download the samples and Prepare sample info.
cd $work_dir
# download dataset from >> https://uploadfiles.io/kc0qqbvd 
unzip ngs2-assignment-data.zip

R1="$work_dir/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq.gz"
R2="$work_dir/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq.gz"

R1s="$work_dir/ngs2-assignment-data/shuffled_SRR8797509_1.part_001.part_001.fastq.gz"
R2s="$work_dir/ngs2-assignment-data/shuffled_SRR8797509_2.part_001.part_001.fastq.gz"

SM="SRR8797509"		# sample name
LB=$SM				# library name
PL="Illumina"		# Platform
PU="HiSeqXTen"	# platform unit
RGID="null"			# Read-Group ID

#------------------------------------------------------------
# 3- Map the sample to the ref. (output=SAM)
mkdir hisatmap 
hisat2 -q --phred33 -x $index -1 $R1 -2 $R2 -S hisatmap/$SM.map.sam --met-file hisatmap/mapping.met > hisatmap/map.log

#------------------------------------------------------------
# 4- Add read groups, sort, mark duplicates, and create index (Picard tool) (input=sam; output=dedup.bam & index)

mkdir dedup

java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=hisatmap/$SM.map.sam O=dedup/$SM.rg_added.sorted.bam SO=coordinate RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM

java -Xmx2g -jar $picard_path/picard.jar MarkDuplicates I=dedup/$SM.rg_added.sorted.bam O=dedup/$SM.dedup.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=dedup/output.met

#------------------------------------------------------------
# 5- Split'N'Trim and reassign mapping qualities
mkdir splitN

gatk SplitNCigarReads -R $ref_fna -I dedup/$SM.dedup.bam -O splitN/$SM.splitN.bam

#------------------------------------------------------------
# 6- Base Recalibration
# m4 ha3rf a3mlo 34an m4 m3aya vcf file bel known SNP.
mkdir BR 
# ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz	-o BR/00-All.vcf.gz	# it's 16.3 GB!
gunzip 00-All.vcf.gz
vcf_file="$work_dir/BR/00-All.vcf"

gatk --java-options "-Xmx2G" BaseRecalibrator -R $ref_fna -I splitN/$SM.splitN.bam --known-sites $vcf_file -O BR/$SM.report

gatk --java-options "-Xmx2G" ApplyBQSR -R $ref_fna -I splitN/$SM.splitN.bam -bqsr BR/$SM.report -O BR/$SM.bqsr.bam --add-output-sam-program-record --emit-original-quals

#------------------------------------------------------------
# 8- Variant calling
mkdir $work_dir/vc
gatk --java-options "-Xmx2G" HaplotypeCaller -R $ref_fna -I splitN/$SM.splitN.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0 -O vc/$SM.vcf
# when BR works, use its output instead!


#------------------------------------------------------------
# 9- Variant filtering
gatk --java-options "-Xmx2G" VariantFiltration -R $ref_fna -V vc/$SM.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O vc/$SM.final.vcf 

