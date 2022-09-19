#!/usr/bin/bash 
# /where_is_tools==''
conda activate WGS
module purge
module load singularity
id=$1
#trim the fastq file
trimmomatic PE 
   /where_is_your_fastq_files/${id}_R1_001.fastq.gz 
   /where_is_your_fastq_files/${id}_R2_001.fastq.gz 
   /where_is_your_data/$id/fq/${id}_R1_001.paired.fastq.gz 
   /where_is_your_data/$id/fq/${id}_R1_001.unpaired.fastq.gz 
   /where_is_your_data/$id/fq/${id}_R2_001.paired.fastq.gz 
   /where_is_your_data/$id/fq/${id}_R2_001.unpaired.fastq.gz 
   ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# bam file
singularity run /where_is_tools/bwa.img bwa mem -t 8 -M -Y 
   -R "@RG\tID:S200\tPL:ILLUMINA\tLB:L008\tSM:${id}" 
   /where_is_your_references/mm10/chr1_19_X_Y_M/GRCm38.chr1_19_X_Y_M.fa 
   /where_is_your_data/$id/fq/${id}_R1_001.paired.fastq.gz 
   /where_is_your_data/$id/fq/${id}_R2_001.paired.fastq.gz 
   | singularity run ../../tools/samtools.img samtools view -Sb 
   - > /where_is_your_data//$id/bam/${id}.bam

# sorted bam
singularity run /where_is_tools/samtools.img samtools sort -@ 4 -m 4G -O bam 
   -o /where_is_your_data/$id/bam/${id}.sorted.bam 
   /where_is_your_data/$id/bam/${id}.bam
sleep 20s
# mark duplicate
singularity run /where_is_tools/gatk.img gatk MarkDuplicates 
    -I /where_is_your_data/${id}/bam/${id}.sorted.bam 
    -M /where_is_your_data/${id}/bam/${id}.markdup_metrics.txt 
    -O /where_is_your_data/${id}/bam/${id}.sorted.markdup.bam
sleep 20s
# index bam
singularity run /where_is_tools/samtools.img samtools index /where_is_your_data/${id}/bam/${id}.sorted.markdup.bam
##BQSR
singularity run /where_is_tools/gatk.img gatk BaseRecalibrator 
     -R /where_is_your_references/mm10/chr1_19_X_Y_M/GRCm38.chr1_19_X_Y_M.fa
     -I /where_is_your_data/${id}/bam/${id}.sorted.markdup.bam 
     --known-sites /where_is_your_references/mm10/GATK_mm10_bundle/mgp.v5.indels.pass.chr.vcf 
     -O /where_is_your_data/${id}/bam/${id}.sorted.markdup.recal_data.table && echo "* ${id}.sorted.markdup.recal_data.table done*"
# Apply the BQSR
singularity run /where_is_tools/gatk.img gatk ApplyBQSR
    --bqsr-recal-file /where_is_your_data/${id}/bam/${id}.sorted.markdup.recal_data.table 
    -R /where_is_your_references/mm10/chr1_19_X_Y_M/GRCm38.chr1_19_X_Y_M.fa
    -I /where_is_your_data/${id}/bam/${id}.sorted.markdup.bam
    -O /where_is_your_data/${id}/bam/${id}.sorted.markdup.BQSR.bam && echo "* ApplyBQSR done*"

#index for ${sample}.sorted.markdup.BQSR.bam
singularity run ../where_is_tools/samtools.img samtools index 
     /where_is_your_data/${id}/bam/${id}.sorted.markdup.BQSR.bam
     && echo "*${sample}.sorted.markdup.BQSR.bam index done *"
#end
