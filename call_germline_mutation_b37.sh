# BAM file used for variant calling
bamfile=$1
# Sample id
sid=$2
# Target regions
target=$3

##################################### Please modify the paths/database files ##############################################
# adjust the following path to reflect where you have downloaded the dgvar package
workdir=/code/dgvar-master

# adjust the following paths to reflect where you installed the tools
# samtools
samtools=/samtools-1.9/samtools
module load singularity
#module load java
# java v1.6 (used for GATK)
java6="singularity run /code/dgvar-master/bjava6.img java"
# java v1.7 (used for snpEff/snpSift)
java7=java
#"singularity run /code/dgvar-master/java7.img java"
# python v2.7+
python=python
# GATK v2.5.2
gatk=/GenomeAnalysisTK/GenomeAnalysisTK.jar

# sneEff v4.2
snpeff=/snpEff/snpEff.jar
snpsift=/snpEff/SnpSift.jar
config=/snpEff/snpEff.config

# databases
dbdir=${workdir}/db
# human hg19 reference genome, please make sure you also have the dictionary file (.dict) and index file (.fai) in the save path
refseq=${dbdir}/human_g1k_v37.fasta
# dbSNP hg19 v137
snpdb=${dbdir}/dbsnp_138.b37.vcf
# ClinVar v20180805 (please download from our dgvar github repository)
clinvar=${dbdir}/clinvar_20180805.fix.vcf.gz
# ExAC database (please download from our dgvar github repository)
exac=${dbdir}/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.PASS.vcf.gz
exacnopass=${dbdir}/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.nonPASS.vcf.gz
# Gene annotations (Onco genes and Tumor Suppressor genes, please download from our dgvar github repository)
geneann=${dbdir}/gene_annotations.txt
# Common variants in the EIPM cohort (wnload from our dgvar github repository)
freqtablefile=${dbdir}/eipm_common_variants.txt
##########################################################################################################################

# configuration
mem=10g
anndb=GRCh37.75
flag="-no-upstream -no-downstream -no-intergenic -v"
ExACcols="AF,AdjAF,AN_Adj,AF_AFR,AF_AMR,AF_EAS,AF_FIN,AF_NFE,AF_SAS,AF_OTH,AF_MALE,AF_FEMALE"
CLNcols="CLNVARID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID,ORIGIN,GENEINFO,MC,RS,SSR"

# data folders
srcdir=${workdir}/src
outdir=${workdir}/results
vardir=${outdir}/${sid}
varlogdir=${vardir}/logs
vartmpdir=${vardir}/var_tmp
anndir=${vardir}/annotated
annlogdir=${anndir}/logs
anntmpdir=${anndir}/ann_tmp
filtdir=${anndir}/filt
filtlogdir=${filtdir}/logs

####################### functions #######################
# check whether the previous command was successful
# if not, exit with a pre-defined error message
mycheck(){
	if [ $? -ne 0 ]
	then
		echo -e "\nError: $1"
		exit 99
	fi
}
#########################################################

echo "Start to process Sample ${sid}."

# check output folder
if [ -d ${vardir} ]
then
	echo "Output folder already exists, nothing is done."
	exit 0
else
	mkdir ${vardir}
	mkdir ${varlogdir}
	mkdir ${vartmpdir}
	mkdir ${anndir}
	mkdir ${annlogdir}
	mkdir ${anntmpdir}
	mkdir ${filtdir}
	mkdir ${filtlogdir}
fi

# check source bam
echo -n "Checking source bam file..."
if [ -e ${bamfile} ]
then
        echo "ok."
else
        echo -e "\nCannot locate the source bam file: ${bamfile}"
	exit 1
fi

# copy source bam (and .bai) file to local directory
echo -n "Copying bam and preparing idx bai..."
cd ${vardir}
mycheck "Failed to change to directory ${vardir}"

bam=Sample_${sid}.bam
cp ${bamfile} ${bam}
mycheck "Failed to copy source bam ${bamfile} to the local directory `pwd`"

bamprefix=${bamfile%.bam}
if [ -e ${bamprefix}.bai ]
then
	cp ${bamprefix}.bai  Sample_${sid}.bai
elif [ -e ${bamfile}.bai ]
then
	cp ${bamfile}.bai  ${bam}.bai
else
	${samtools} index ${bam}
fi
mycheck "Failed to get index for bam ${bamfile}"
echo "ok."

# produce raw variant calls
echo -n "Make raw variant calls using GATK UnifiedGenotyper..."
${java6} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -glm BOTH -R ${refseq} -T UnifiedGenotyper -D ${snpdb} -o ${vardir}/${sid}.UG.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -minIndelFrac 0.1 -dcov 5000 -A AlleleBalance -A BaseCounts -A GCContent -A LowMQ -A SampleList -A VariantType -L ${target} -nt 8 -nct 4 \
-I ${bam} \
>${varlogdir}/${sid}.GATK.UnifiedGenotyper.log 2>&1
mycheck "Failed to make variant calls using GATK UnifiedGenotyper"
echo "ok."
# filter raw variant calls
echo -n "Filter raw variants using GATK VariantFiltration..."
${java6} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -R ${refseq} -T VariantFiltration -o ${vardir}/${sid}.UG.filtered.vcf --variant ${vardir}/${sid}.UG.raw.vcf --clusterWindowSize 10 --clusterSize 3  --filterExpression "(DP - MQ0) < 10" --filterName "LowCoverage" \
>${varlogdir}/${sid}.GATK.VariantFiltration.log 2>&1
mycheck "Failed to pre-filter variant calls using GATK VariantFiltration"
echo "ok."

#echo -n "Make raw variant calls using GATK HaplotyperCaller..."
#${gatk} --java-options "-Xmx${mem}" HaplotypeCaller -R ${refseq} -D ${snpdb} -O ${vardir}/${sid}.UG.raw.vcf --standard-min-confidence-threshold-for-calling 50.0 -L ${target}  \
#-I ${bam} \
#>${varlogdir}/${sid}.GATK.HaplotypeCaller.log 2>&1
#mycheck "Failed to make variant calls using GATK HaplotyperCaller"
#echo "ok."

# filter raw variant calls
#echo -n "Filter raw variants using GATK VariantFiltration..."
#${gatk} --java-options "-Xmx${mem}" VariantFiltration -R ${refseq} -O ${vardir}/${sid}.UG.filtered.vcf --variant ${vardir}/${sid}.UG.raw.vcf --cluster-window-size 10 --cluster-size 3  --filter-expression "(DP - MQ0) < 10" --filter-name "LowCoverage" \
#>${varlogdir}/${sid}.GATK.VariantFiltration.log 2>&1
#mycheck "Failed to pre-filter variant calls using GATK VariantFiltration"
#echo "ok."

# annotate variants using snpeff
echo -n "Annotate variants using snpEff..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpeff} ann -c ${config} -s ${anndir}/${sid}.summary.snpeff.UG.html ${flag} ${anndb} ${vardir}/${sid}.UG.filtered.vcf >${anndir}/${sid}.UG.filtered.snpeff.vcf 2>${annlogdir}/${sid}.snpeff.ann.log
mycheck "Failed to annotate variants using snpEff"
echo "ok."

# add clinVar annotations using snpsift
echo -n "Add ClinVar annotations using snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name clinvar -info ${CLNcols} ${clinvar} ${anndir}/${sid}.UG.filtered.snpeff.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.vcf 2>${annlogdir}/${sid}.snpsift.clinvar.log
mycheck "Failed to add ClinVar annotations"
echo "ok."

# add ExAC AF info using snpsift
echo -n "Add ExAC AF (allele frequency) annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exac} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.log
mycheck "Failed to add ExAC AF annotations"
#echo "ok."

# add ExAC filterings using snpsift
echo -n "Add ExAC filterings by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.extend.log
mycheck "Failed to add ExAC filterings"
echo "ok."
gzip ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf
mycheck "Failed to gzip vcf ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf"

# update FILTER column
echo -n "Set ExAC filtering..."
${python} ${srcdir}/set_ExAC_filter.py ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf.gz ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.patch.vcf.gz >${annlogdir}/${sid}.set_ExAC_filter.log
mycheck "Failed to set ExAC filtering in vcf"
echo "ok."

# extract info and screen
echo -n "Pre-categorize variants..."
${python} ${srcdir}/screen_germline_variants.py ${sid} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.patch.vcf.gz ${filtdir}/ ${filtdir}/${sid}.screen_variants.warning.txt >${filtlogdir}/${sid}.screen_germline_variants.log 2>&1
mycheck "Failed to pre-categorize variants"
echo "ok."

# add gene annotations
echo -n "Add gene annotation..."
${python} ${srcdir}/add_gene_annotations.py ${filtdir}/${sid}.info.all.txt.gz ${geneann} ${filtdir}/${sid}.info.all.ann.txt.gz >${filtlogdir}/${sid}.add_gene_annotations.log 2>&1
mycheck "Failed to add gene annotation"
echo "ok."

# select candidates
echo -n "Selecting candidate variants..."
${python} ${srcdir}/select_candidate_variants.py ${sid} ${filtdir}/${sid}.info.all.ann.txt.gz ${filtdir}/${sid}.candidates.txt.gz >${filtlogdir}/${sid}.select_candidate_variants.log 2>&1
mycheck "Failed to select candidate variants"
echo "ok."

# filter common variants in the EIPM cohort
echo -n "Filter common variants..."
${python} ${srcdir}/filter_common_variants.py ${filtdir}/${sid}.candidates.txt.gz ${freqtablefile} ${filtdir}/${sid}.candidates.filter_common.txt.gz >${filtlogdir}/${sid}.filter_common_variants.log 2>&1
mycheck "Failed to filter common variants"
echo "ok."

# counting candidate variants
numLines=`zcat ${filtdir}/${sid}.candidates.filter_common.txt.gz |wc -l`
numVars=`expr ${numLines} - 1`
echo "${numVars} candidate germline variants detected."

# cleaning
rm -rf ${vartmpdir}
rm -rf ${anntmpdir}
rm -f ${vardir}/${bam}*
rm ${vardir}/${sid}.UG.raw.vcf*
gzip ${vardir}/${sid}.UG.filtered.vcf*
rm -f ${anndir}/${sid}.UG.filtered.snpeff.vcf
rm -f ${anndir}/${sid}.UG.filtered.snpeff.clinvar.vcf
rm -f ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.vcf
rm -f ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf.gz
rm -f ${anndir}/${sid}.summary.snpeff.UG.*
rm -f ${filtdir}/${sid}.info.all*.txt.gz
rm -f ${filtdir}/${sid}.candidates.txt.gz

echo "Complete processing Sample ${sid}."
