#!/bin/bash
module load java
trimmomatic=/path_to_tools/Trimmomatic-0.38/trimmomatic-0.38.jar
bwa=/path_to_tools/bwa-0.7.17/bwa
samtools=/path_to_tools/samtools-1.9/samtools
gatk=/path_to_tools/gatk-4.1.4.0/gatk

#reference
reference=/path_to_reference/b37
GATK_bundle=/path_to_gatk_bundle/GTAK_bundle

#make sure tbi index for bundle
#shell parameters
sample=$1 ## sample id

fq1=$(ls /scratch/fs446/hWGS_ptc/data/$sample/*.1.fq.gz)
fq2=$(ls /scratch/fs446/hWGS_ptc/data/$sample/*.2.fq.gz)
##make fastq1 and fastq2 files are the same prefix, such as *.1.fq.gz and *.2.fq.gz

fq_file_name=`basename $fq1`
fq_file_name=${fq_file_name%%.1.fq.gz}
##eg. fq_file_name=TP002_CSFP190048404-1a_HTVFWDSXX_L4

RGID="${fq_file_name:23:12}" ## lan id
library="${fq_file_name:6:16}" ## library id
outdir=/scratch/fs446/hWGS_ptc/output

##set folder by sample id
outdir=${outdir}/${sample}

#output directory

if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

#if [ ! -d $outdir/gatk ]
#then mkdir -p $outdir/gatk
#fi

##using Trimmomatic for data QC, set keepBothReads = True in ILLUMINACLIP
time java -jar ${trimmomatic} PE 
  $fq1 
  $fq2 
  $outdir/cleanfq/${fq_file_name}.paired.1.fq.gz ${fq_file_name}.unpaired.1.fq.gz 
  $outdir/cleanfq/${fq_file_name}.paired.2.fq.gz ${fq_file_name}.unpaired.2.fq.gz 
  ILLUMINACLIP:/home/fs446/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:True SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 &&
  echo "*fq QC done *"

##using bwa mem for alignment
time $bwa mem -t 8 -M -Y 
  -R "@RG\tID:$RGID\tPL:ILLUMINA\tLB:$library\tSM:$sample" 
  $reference/human_g1k_v37.fasta $outdir/cleanfq/${fq_file_name}.paired.1.fq.gz 
  $outdir/cleanfq/${fq_file_name}.paired.2.fq.gz | $samtools view -Sb - > $outdir/bwa/${sample}.bam && 
  echo "* BWA MEM done *" && time $samtools sort -@ 4 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "* sorted raw bam file done * "

#mark duplicates
time $gatk MarkDuplicates 
  -I $outdir/bwa/${sample}.sorted.bam 
  -M $outdir/bwa/${sample}.markdup_metrics.txt 
  -O $outdir/bwa/${sample}.sorted.markdup.bam && 
  echo "* ${sample}.sorted.bam MarkDuplicates done *"

## index for ${sample}.sorted.markdup.bam
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && 
  echo " *${sample}.sorted.markdup.bam index done *"

##BQSR
# make sure references have indexed

time $gatk BaseRecalibrator -R $reference/human_g1k_v37.fasta -I $outdir/bwa/${sample}.sorted.markdup.bam 
  --known-sites $GATK_bundle/1000G_phase1.snps.high_confidence.b37.vcf 
  --known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf 
  --known-sites $GATK_bundle/dbsnp_138.b37.vcf 
  -O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && 
  echo "* ${sample}.sorted.markdup.recal_data.table done*"

time $gatk ApplyBQSR 
  --bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table 
  -R $reference/human_g1k_v37.fasta 
  -I $outdir/bwa/${sample}.sorted.markdup.bam 
  -O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && 
  echo "* ApplyBQSR done*"

#index for ${sample}.sorted.markdup.BQSR.bam
time $samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && 
echo "*${sample}.sorted.markdup.BQSR.bam index done *"

$gatk HaplotypeCaller 
  --emit-ref-confidence GVCF 
  -R $reference/human_g1k_v37.fasta 
  -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam 
  -O $outdir/gatk/${sample}.HC.g.vcf.gz && 
  echo "* GVCF ${sample}.HC.g.vcf.gz done *"

### SNP mode
time $gatk VariantRecalibrator 
  -R $reference/human_g1k_v37.fasta 
  -V $outdir/population/${outname}.HC.vcf.gz 
  -AS -resource:hapmap,known=false,training=true,truth=true,prior=15.0 
  $GATK_bundle/hapmap_3.3.b37.vcf 
  -resource:omini,known=false,training=true,truth=false,prior=12.0 
  $GATK_bundle/1000G_omni2.5.b37.vcf 
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 
  $GATK_bundle/1000G_phase1.snps.high_confidence.b37.vcf 
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 
  $GATK_bundle/dbsnp_138.b37.vcf 
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP 
  -tranche 100.0 -tranche 99.0 -tranche 95.0 -tranche 92.0 -tranche 90.0  
  --max-gaussians 6 -rscript-file $outdir/population/${outname}.HC.snps.plot.R 
  --tranches-file $outdir/population/${outname}.HC.snps.tranches 
  -O $outdir/population/${outname}.HC.snps.recal

sleep 100s

time $gatk ApplyVQSR 
  -R $reference/human_g1k_v37.fasta 
  -V $outdir/population/${outname}.HC.vcf.gz 
  -ts-filter-level 99.0 
  --tranches-file $outdir/population/${outname}.HC.snps.tranches 
  -recal-file $outdir/population/${outname}.HC.snps.recal 
  -mode SNP 
  -O $outdir/population/${outname}.HC.snps.VQSR.vcf.gz && 
  echo "* SNPs VQSR done*"

# INDEL mode
time $gatk VariantRecalibrator 
  -R $reference/human_g1k_v37.fasta 
  -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz 
  -resource:mills,known=true,training=true,truth=true,prior=12.0 
  $GATK_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf 
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL 
  --max-gaussians 6 
  -rscript-file $outdir/population/${outname}.HC.snps.indels.plot.R 
  --tranches-file $outdir/population/${outname}.HC.snps.indels.tranches 
  -O $outdir/population/${outname}.HC.snps.indels.recal && 
  echo "*recal finished"

sleep 100s
time $gatk ApplyVQSR 
  -R $reference/human_g1k_v37.fasta 
  -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz 
  -ts-filter-level 99.0 
  --tranches-file $outdir/population/${outname}.HC.snps.indels.tranches 
  -recal-file $outdir/population/${outname}.HC.snps.indels.recal 
  -mode INDEL 
  -O $outdir/population/${outname}.HC.VQSR.vcf.gz && 
  echo "* SNPs and Indels VQSR (${outname}.HC.VQSR.vcf.gz) finsih) done****"
