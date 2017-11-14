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