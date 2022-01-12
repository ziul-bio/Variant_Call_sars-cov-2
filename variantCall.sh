#!/usr/bin/env bash

# Para terminar a execução do script se algum erro ocorrer
set -o errexit

conda activate genomics

# ------------------------------------------------------------------------------------------------------------------------------

#                                         VARIANT CALLING

#-------------------------------------------------------------------------------------------------------------------------------                               

# Esta atividade simula a execução de um pipeline simplificado para identificação de variantes (genotipagem) do vírus SARS-CoV-2



# ---------------------------- Organizando o diretório de trabalho - pwd ----------------------------

# Criando os diretórios, onde os arquivos serão armazenados
mkdir -p ref_genome raw_fastq passedQC mapped variants


# ---------------------------- Download genoma de referência e anotação ------------------------------

# Genoma de referência e anotação
wget -P ref_genome/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
wget -P ref_genome/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz

# Para extrair os arquivos use:
gzip -d ref_genome/GCF_009858895.2_ASM985889v3_genomic.fna.gz
gzip -d ref_genome/GCF_009858895.2_ASM985889v3_genomic.gff.gz

# Para renomear use:
mv ref_genome/GCF_009858895.2_ASM985889v3_genomic.fna ref_genome/covid_genome.fna
mv ref_genome/GCF_009858895.2_ASM985889v3_genomic.gff ref_genome/covid_genome.gff


# Download dos arquivos fastq provenientes do sequenciamento
wget -P raw_fastq/ http://projetos.lbi.iq.usp.br/phaghost/data/FASTQ/covid_patient_R1.fastq.gz
wget -P raw_fastq/ http://projetos.lbi.iq.usp.br/phaghost/data/FASTQ/covid_patient_R2.fastq.gz



# -------------------------- Controle de Qualidade das Sequências com Trimomatic ---------------------

# Checando a qualidade das reads antes da trimagem
fastqc raw_fastq/*.fastq.gz

# Trimmomatic trimagem das reads

fastp -q 30 -l 40 -f 6 -F 6 -t 9 -T 9 -i raw_fastq/covid_patient_R1.fastq.gz -I raw_fastq/covid_patient_R2.fastq.gz -o passedQC/covid_patient_R1.fastq.gz -O passedQC/covid_patient_R2.fastq.gz

# Checando as reads pós CQ.
fastqc passedQC/*.fastq.gz


# ------------------ Mapeamento dos reads no genoma de referência do SARS-CoV-2 ----------------------

# Gerando index do arquivo .fa > .fai
bwa index -p ref_genome/covid ref_genome/covid_genome.fna

# Alinhamento, ordenamento por nomes
bwa mem -M -t 4 ref_genome/covid passedQC/covid_patient_R1.fastq.gz passedQC/covid_patient_R2.fastq.gz | samtools sort -n - > mapped/covid_patient.sorted_Names.bam 

# This command takes the file sorted by name resulted from alignment and run in pipe samtools fixmate, samtools sort(now by coordenate) and samtools markdup
samtools fixmate -m mapped/covid_patient.sorted_Names.bam - | samtools sort -o - | samtools markdup -r -s - mapped/covid_patient.sorted.Markdup.bam

# Criando um arquivo para checar estatísticas do alinhamento com samtools flagstat
samtools flagstat mapped/covid_patient.sorted.Markdup.bam



# ---------------------------- Variant Calling with Freebayes ------------------------------------------

# Gerando index do arquivo .bam > .bai
samtools index mapped/covid_patient.sorted.Markdup.bam

# Executando Freebayes e Filtrando o arquivo vcf
freebayes -p2 -f ref_genome/covid_genome.fna mapped/covid_patient.sorted.Markdup.bam | bcftools filter -e "QUAL<20 || GT!='1/1'" - -o variants/covid_patient.filtered.vcf


# Variant Annotation
snpEff -v -classic NC_045512.2 variants/covid_patient.filtered.vcf > variants/covid_patient.snpEff.vcf


# Subseting the filds of interest
cat variants/covid_patient.snpEff.vcf | cut -f8 | cut -d "|" -f4,6 > variants/AA_chane.txt



