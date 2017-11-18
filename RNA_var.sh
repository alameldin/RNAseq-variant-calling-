## Define the path of your raw data and your scripts
data_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR/data"
script_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR/scripts"
out_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR/output"
genome_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR/genomeDir"
#####################################
##Downloading the genome data:

cd ${genome_path}
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r77.all.fa
# downloading the annotation file
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip Homo_sapiens.GRCh38.77.gtf.gz

####
#running the STAR_index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${genome_path}/GRCh38 star_indices_overhang100/ --genomeFastaFiles ${genome_path}/GRCh38_r77.all.fa --sjdbGTFfile ${genome_path}/Homo_sapiens.GRCh38.77.gtf --sjdbOverhang 100
#########
# running the star alignement 

STAR --genomeDir /mnt/ls15/scratch/users/hussien/RNA_VAR/genomeDir/GRCh38 --readFilesIn ${data_path}/ERR1050075_1.fastq ${data_path}/ERR1050075_2.fastq  --runThreadN 7

#######
# For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass

STAR --runMode genomeGenerate --genomeDir ${genome_path}/new_index --genomeFastaFiles ${genome_path}/GRCh38_r77.all.fa  --sjdbFileChrStartEnd ${out_path}/star_out/SJ.out.tab --sjdbOverhang 100 --runThreadN 4
##########

module load Java/1.8.0_31
module load  picardTools/1.89
java -jar $PICARD/AddOrReplaceReadGroups.jar I=/mnt/ls15/scratch/users/hussien/RNA_VAR/output/new_star_out/Aligned.out.sam O=${genome_path}/new_index/rg_added_sorted.bam Scoordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
java -jar $PICARD/MarkDuplicates.jar I=${genome_path}/new_index/rg_added_sorted.bam O=${genome_path}/new_index/dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

######
#Split'N'Trim and reassign mapping qualities
# we need to prepare the .fasta genome file to be used as a reference in GATK first using either SAMtools and PICARDtools.
mkdir ${genome_path}/GATK_indexed
cd ${genome_path}/GATK_indexed
module load SAMtools
samtools faidx ${genome_path}/GRCh38_r77.all.fa

### creathing .dict file for the reference genome unsing Picardtools
java -jar $PICARD/CreateSequenceDictionary.jar R= ${genome_path}/GRCh38_r77.all.fa O= ${genome_path}/GATK_indexed/GRCh38_r77.all.dict
#######
#Split'N'Trim and reassign mapping qualities
module load GATK
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -I ${genome_path}/new_index/dedupped.bam -o ${out_path}/GATK_out/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#####
#Variant calling

java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -I ${out_path}/GATK_out/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${out_path}/GATK_out/output.vcf

######
# Variant filtering
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -V ${out_path}/GATK_out/output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${out_path}/GATK_out/output_filtered.vcf

