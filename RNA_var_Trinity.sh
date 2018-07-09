## Define the path of your raw data and your scripts
data_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_Trinity/data"
script_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_Trinity/scripts"
out_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_Trinity/output"
genome_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_Trinity/genomeDir"
#####################################
##Downloading the RNA data:
cd ${data_path}
module load GNU/4.4.5
module load SRAToolkit/2.8.2
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1258218
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1153470

########
cd ${genome_path}
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r77.all.fa
# downloading the annotation file
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip Homo_sapiens.GRCh38.77.gtf.gz
###################
##Running the Trinity
cd ${script_path}
qsub Trinity_470.sh
qsub Trinity_218.sh
####
#running the STAR_index
mkdir GRCh38
qsub ${script_path}/star_index.sh
#STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${genome_path}/GRCh38 star_indices_overhang100/ --genomeFastaFiles ${genome_path}/GRCh38_r77.all.fa --sjdbGTFfile ${genome_path}/Homo_sapiens.GRCh38.77.gtf --sjdbOverhang 100
#########
# running the star alignement 
mkdir ${genome_path}/STAR_mapping
qsub ${script_path}/star_mapping.sh
#STAR --genomeDir ${genome_path}/GRCh38 --readFilesIn ${data_path}/SRR1153470_1.fastq ${data_path}/SRR1153470_2.fastq  --runThreadN 7
#######

# For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
qsub ${script_path}/star_index_1.sh
#STAR --runMode genomeGenerate --genomeDir ${genome_path}/new_index --genomeFastaFiles ${genome_path}/GRCh38_r77.all.fa --sjdbFileChrStartEnd ${genome_path}/STAR_mapping/SJ.out.tab --sjdbOverhang 75 --runThreadN 4
mkdir ${genome_path}/STAR_mapping_1
qsub star_mapping_1.sh

##########
# Mark the duplicates
mkdir ${genome_path}/Picard_index
module load Java/1.8.0_31
module load  picardTools/1.89
java -jar $PICARD/AddOrReplaceReadGroups.jar I=${genome_path}/STAR_mapping_1/Aligned.out.sam O=${genome_path}/Picard_index/rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
java -jar $PICARD/MarkDuplicates.jar I=${genome_path}/Picard_index/rg_added_sorted.bam O=${genome_path}/Picard_index/dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

######
#Split'N'Trim and reassign mapping qualities
# we need to prepare the .fasta genome file to be used as a reference in GATK first using either SAMtools and PICARDtools.
mkdir ${genome_path}/GATK_indexed
cd ${genome_path}/GATK_indexed
cp ${genome_path}/GRCh38_r77.all.fa .
module load SAMtools
samtools faidx GRCh38_r77.all.fa > GRCh38_r77.all.fa.fai

### creathing .dict file for the reference genome unsing Picardtools
java -jar $PICARD/CreateSequenceDictionary.jar R= ${genome_path}/GRCh38_r77.all.fa O= ${genome_path}/GATK_indexed/GRCh38_r77.all.dict
#######
#Split'N'Trim and reassign mapping qualities
mkdir ${out_path}/GATK_out
module load GATK
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -I ${genome_path}/Picard_index/dedupped.bam -o ${out_path}/GATK_out/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
#####Creating Recalibrated BAM file
#java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -I ${out_path}/GATK_out/split.bam -BQSR ${out_path}/GATK_out/recalibration_report.grp -o ${out_path}/GATK_out/recalibrated.bam
#####
#Variant calling

java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -I ${out_path}/GATK_out/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${out_path}/GATK_out/output.vcf

########
module load BEDTools/2.24.0
samtools view -b ${out_path}/GATK_out/bamout.bam | genomeCoverageBed -ibam  stdin -g ${genome_path}/GRCh38_r77.all.fa -bg -split > ${out_path}/GATK_out/bed_out_split.bed

sed 's/^/chr/' ${out_path}/GATK_out/bed_out_split.bed > ${out_path}/GATK_out/bed_out_split_CorID.bed
##### This step was to filter the bed file based on coverage (to elemenate the variants that were in areas covered less than 1 or 5)
cat ${out_path}/GATK_out/bed_out_split_CorID.bed | awk '$4>=1 {print $1 "\t" $2 "\t" $3}' > ${out_path}/GATK_out/out_split_CorID_1.bed
cat ${out_path}/GATK_out/bed_out_split_CorID.bed | awk '$4>=5 {print $1 "\t" $2 "\t" $3}' > ${out_path}/GATK_out/out_split_CorID_5.bed

#to cluster and categorize the bed file
bedtools cluster -i ${out_path}/GATK_out/out_split_CorID_1.bed > ${out_path}/GATK_out/out_split_CorID_clustered_1.bed
bedtools cluster -i ${out_path}/GATK_out/out_split_CorID_5.bed > ${out_path}/GATK_out/out_split_CorID_clustered_5.bed

# bed files manipulations
#multiIntersectBed -i /mnt/ls15/scratch/users/hussien/RNA_VAR_NEW/Genome_var/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed out_split_CorID_clustered_1.bed > compair_GvsR_1.bed


#vcftools --gzvcf myfile.vcf.gz --recode --bed bedfile.bed --out outfile
module load vcftools

#sed 's/^/chr/' ${out_path}/GATK_out/out_split_clustered_1.bed > ${out_path}/GATK_out/out_split_clustered_1_CorID.bed
#sed 's/^/chr/' ${out_path}/GATK_out/out_split_clustered_5.bed > ${out_path}/GATK_out/out_split_clustered_5_CorID.bed
vcftools --gzvcf  ${genome_var}/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz --recode --bed ${out_path}/GATK_out/out_split_CorID_clustered_1.bed --out ${out_path}/GATK_out/out_vcf_1.vcf

vcftools --gzvcf  ${genome_var}/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz --recode --bed ${out_path}/GATK_out/out_split_CorID_clustered_5.bed --out ${out_path}/GATK_out/out_vcf_5.vcf
######
# Variant filtering
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -V ${out_path}/GATK_out/output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${out_path}/GATK_out/output_filtered.vcf

########
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz.tbi

#########
Comparing the RNA variant calling output with the already published DNA variant calling output (downloaded last 2 lines)
mkdir /mnt/ls15/scratch/users/hussien/RNA_VAR/vcftool
cd /mnt/ls15/scratch/users/hussien/RNA_VAR/vcftool
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz.tbi
mv NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz Genome_var.vcf
cp ${out_path}/output.vcf RNA_var.vcf

###
#indexing the vcf file
module load tabix
module load vcftools

bgzip -c RNA_var.vcf > RNA_var.vcf.gz
bgzip -c Genome_var.vcf > Genome_var.vcf.gz
tabix -p vcf *.gz
####
# some statistics on the comparizon between the two files
vcf-compare  Genome_var.vcf.gz RNA_var.vcf.gz > compairing_output_4

more compairing_output_4
# This file was generated by vcf-compare.
# The command line was: vcf-compare(r940) Genome_var.vcf.gz RNA_var.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are:
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
#VN    2221    Genome_var.vcf.gz (0.0%)    RNA_var.vcf.gz (1.4%)
#VN    151550    RNA_var.vcf.gz (98.6%)
#VN    14559593    Genome_var.vcf.gz (100.0%)
###SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
#SN    Number of REF matches:    1513
#SN    Number of ALT matches:    1393
#SN    Number of REF mismatches:    708
#SN    Number of ALT mismatches:    56
#SN    Number of samples in GT comparison:    0
# Number of sites lost due to grouping (e.g. duplicate sites): lost, %lost, read, reported, file
SN    Number of lost sites:    7623    0.1%    14569437    14561814    Genome_var.vcf.gz


#####
#split the variant that is in the genome but not in the transcriptome
vcf-isec -c -f Genome_var.vcf.gz RNA_var.vcf.gz > compairing.vcf_2

#split the variant that is in the Transcriptome but not in the Genome
vcf-isec -c -f RNA_var.vcf.gz Genome_var.vcf.gz > compairing.vcf_3

#####
it looks like the interaction is only in 2197 positions
wc -l compairing.vcf_2
# 14567358 compairing.vcf_2
wc -l compairing.vcf_3
# 151692 compairing.vcf_3
wc -l NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf
# 14569555 NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf
wc -l RNA_var.vcf
# 153824 RNA_var.vcf


