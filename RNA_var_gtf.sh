## Define the path of your raw data and your scripts
data_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/data"
script_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/scripts"
out_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/output"
genome_path=$"/mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/genomeDir"
#####################################
##Downloading the genome data:

cd ${genome_path}
wget ftp://ftp.completegenomics.com/ReferenceFiles/build37.fa.bz2
bzip2 -d build37.fa.bz2
#This is for generating new assembly file with correct header names
sed 's/^>chr/>/' build37.fa > build37_corID.fa
# downloading the annotation file
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

####
#running the STAR_index
cd ${script_path}
qsub star_index.sh
#########
# running the star alignement 

cd ${script_path}
mkdir ${genome_path}/STAR_mapping
qsub star_mapping.sh
#######
# For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass

cd ${script_path}
qsub star_index_1.sh
##########
# running the star alignement

cd ${script_path}
qsub star_mapping_1.sh
##############################

module load Java/1.8.0_31
module load  picardTools/1.89
mkdir ${genome_path}/Picard_index
java -jar $PICARD/AddOrReplaceReadGroups.jar I=${genome_path}/STAR_mapping_1/Aligned.out.sam O=${genome_path}/Picard_index/Build_37_corID.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
java -jar $PICARD/MarkDuplicates.jar I=${genome_path}/Picard_index/Build_37_corID.bam O=${genome_path}/Picard_index/Build_37_corID_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

######
#Split'N'Trim and reassign mapping qualities
# we need to prepare the .fasta genome file to be used as a reference in GATK first using either SAMtools and PICARDtools.
mkdir ${genome_path}/GATK_indexed
cd ${genome_path}/GATK_indexed
module load SAMtools
samtools faidx ${genome_path}/build37_corID.fa > build37_corID.fa.fai

### creathing .dict file for the reference genome unsing Picardtools
java -jar $PICARD/CreateSequenceDictionary.jar R= ${genome_path}/build37_corID.fa O= ${genome_path}/GATK_indexed/build37_corID.dict
#####
#Variant calling
cd ${genome_path}/GATK_indexed
cp ${genome_path}/build37_corID.fa .
cp ${genome_path}/build37_corID.fa.fai .
mkdir ${out_path}/GATK_out
# Split'N'Trim and reassign mapping qualities
qsub ${script_path}/GATK_1.sh
# Variant Calling
qsub ${script_path}/GATK_2.sh
## Forcing an output in a region that is not covered in the bamout (unnecessarily step only for visualization)
qsub ${script_path}/GATK_forcedbam.sh
######
# Variant filtering
#java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${genome_path}/GATK_indexed/GRCh38_r77.all.fa -V ${out_path}/GATK_out/output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${out_path}/GATK_out/output_filtered.vcf

######
module load BEDTools/2.24.0
samtools view -b ${out_path}/GATK_out/bamout.bam | genomeCoverageBed -ibam  stdin -g ${genome_path}/GRCh38_r77.all.fa -bg -split > ${out_path}/GATK_out/bed_out_split.bed
##### This step was to filter the bed file based on coverage (to elemenate the variants that were in areas covered less than 1 or 5)
cat ${out_path}/GATK_out/bed_out_split.bed | awk '$4>=5 {print $1 "\t" $2 "\t" $3}' > ${out_path}/GATK_out/out_split_1.bed
cat ${out_path}/GATK_out/bed_out_split.bed | awk '$4>=1 {print $1 "\t" $2 "\t" $3}' > ${out_path}/GATK_out/out_split_5.bed

bedtools cluster -i ${out_path}/GATK_out/out_split_1.bed > ${out_path}/GATK_out/out_split_clustered_1.bed
bedtools cluster -i ${out_path}/GATK_out/out_split_5.bed > ${out_path}/GATK_out/out_split_clustered_5.bed

#vcftools --gzvcf myfile.vcf.gz --recode --bed bedfile.bed --out outfile
module load vcftools

vcftools --gzvcf  ${out_path}/Ge/Genome_var.vcf.gz --recode --bed ${out_path}/GATK_out/out_split_clustered_1.bed --out ${out_path}/GATK_out/out_vcf_1.vcf

vcftools --gzvcf  ${out_path}/GATK_out/Genome_var.vcf.gz --recode --bed ${out_path}/GATK_out/out_split_clustered_5.bed --out ${out_path}/GATK_out/out_vcf_5.vcf

#########
#Comparing the RNA variant calling output with the already published DNA variant calling output (downloaded last 2 lines)
mkdir /mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/vcftool
cd /mnt/ls15/scratch/users/hussien/RNA_VAR_gtf/vcftool
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/cg_data/NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz.tbi
mv NA19240_lcl_SRR832874.wgs.COMPLETE_GENOMICS.20130401.snps_indels_svs_meis.high_coverage.genotypes.vcf.gz Genome_var.vcf.gz
gunzip Genome_var.vcf.gz
cp ${out_path}/GATK_out/output.vcf RNA_var.vcf
###
#indexing the vcf file
module load tabix
module load vcftools

bgzip -c RNA_var.vcf > RNA_var.vcf.gz
bgzip -c Genome_var.vcf > Genome_var.vcf.gz
tabix -p vcf Genome_var.vcf.gz
tabix -p vcf RNA_var.vcf.gz
####
# some statistics on the comparizon between the two files
vcf-compare  Genome_var.vcf.gz RNA_var.vcf.gz > comparing_output
more comparing_output
# This file was generated by vcf-compare.
# The command line was: vcf-compare(r940) Genome_var.vcf.gz RNA_var.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are:
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN    35875    RNA_var.vcf.gz (23.0%)
VN    120054    Genome_var.vcf.gz (0.8%)    RNA_var.vcf.gz (77.0%)
VN    14441760    Genome_var.vcf.gz (99.2%)
#SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN    Number of REF matches:    118126
SN    Number of ALT matches:    115843
SN    Number of REF mismatches:    1928
SN    Number of ALT mismatches:    345
SN    Number of samples in GT comparison:    0
# Number of sites lost due to grouping (e.g. duplicate sites): lost, %lost, read, reported, file
SN    Number of lost sites:    7623    0.1%    14569437    14561814    Genome_var.vcf.gz


#####
#split the variant that is in the genome but not in the transcriptome
vcf-isec -c -f Genome_var.vcf.gz RNA_var.vcf.gz > compairing.vcf_2

#split the variant that is in the Transcriptome but not in the Genome
vcf-isec -c -f RNA_var.vcf.gz Genome_var.vcf.gz > compairing.vcf_3

#####
it looks like the interaction is only in a few (~1000) positions
wc compairing.vcf_2
#14449464  144494184 1894840514 compairing.vcf_2
wc -l compairing.vcf_3
#36016 compairing.vcf_3
wc -l Genome_var.vcf
#14569555 Genome_var.vcf
wc -l RNA_var.vcf
# 155982 RNA_var.vcf

# generate the Bam file from the output.vcf

module load BEDTools/2.24.0
bedToBam -ubam -g ${genome_path}/build37_corID.fa.fai -i ${out_path}/GATK_out/output.vcf > ${out_path}/GATK_out/output.bam



