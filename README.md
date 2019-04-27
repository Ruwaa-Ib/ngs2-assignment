# ngs2-assignment

0. Create workdir
```
cd ~/workdir
git clone https://github.com/Ruwaa-Ib/ngs2-assignment.git
cd ngs2-assignment
work_dir="$(pwd)"
```

1. Dataset downloaded from >> https://uploadfiles.io/kc0qqbvd 
```
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
```

2. Download and Index the refrence human genome 
```
# download ref.
mkdir genome-data && cd genome-data
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
java -Xmx2g -jar $picard_path/picard.jar \
  CreateSequenceDictionary \
  R=$ref_fna \
  O=$ref_fna.dict
```

3. Map the dataset to the indexed genome
```
mkdir hisatmap 
hisat2 -q --phred33 \
  -x $index \
  -1 $R1 \
  -2 $R2 \
  -S hisatmap/$SM.map.sam \
  --met-file hisatmap/mapping.met \
  > hisatmap/map.log
```

4. Deduplicate and create index
```
mkdir dedup

java -Xmx2g -jar $picard_path/picard.jar \
  AddOrReplaceReadGroups \
  I=hisatmap/$SM.map.sam \
  O=dedup/$SM.rg_added.sorted.bam \
  SO=coordinate \
  RGID=$RGID \
  RGLB=$LB \
  RGPL=$PL \
  RGPU=$PU \
  RGSM=$SM

java -Xmx2g -jar $picard_path/picard.jar \
  MarkDuplicates \
  I=dedup/$SM.rg_added.sorted.bam \
  O=dedup/$SM.dedup.bam  \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT \
  M=dedup/output.met
```

5. Split'N'Trim
```
mkdir splitN
gatk SplitNCigarReads -R $ref_fna -I dedup/$SM.dedup.bam -O splitN/$SM.splitN.bam
```

6. Base Recalibration

  This part wasn't done becuse the VCF file of the human whole genome is too big to be downloaded and used; plus, its effect is mariginal (according to GATK discussion post).
  
7. Variant Calling

  Somethins wrong is going here! IDK why it produces empty vcf file.
```
mkdir $work_dir/vc
gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R $ref_fna \
  -I splitN/$SM.splitN.bam \
  --dont-use-soft-clipped-bases \
  -stand-call-conf 20.0 \
  -O vc/$SM.vcf
```

8. Variant Filtering
```
gatk --java-options "-Xmx2G" VariantFiltration \
  -R $ref_fna -V vc/$SM.vcf \
  -window 35 \
  -cluster 3 \
  --filter-name FS -filter "FS > 30.0" \
  --filter-name QD -filter "QD < 2.0" \
  -O vc/$SM.final.vcf 
```
