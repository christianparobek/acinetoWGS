### This is a BWA-MEM alignment pipeline
### Modified for Acinetobacter baumannii whole-genome sequencing
### Started December, 2014
### Features:
###	Paired-end alignments using BWA-MEM    
###	Variant-calling using GATK


##########################################################################
###################### REFERENCE GENOME PREPARATION ######################
##########################################################################

## INDEX REFERENCE SEQUENCE FOR BWA (THE ALIGNER)
#bwa index -a bwtsw polished_assembly.fasta

## INDEX SEQUENCE DICTIONARIES FOR GATK (THE REALIGNER / VARIANT CALLER)
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/CreateSequenceDictionary.jar R=polished_assembly.fasta O=polished_assembly.dict

## INDEX REFERENCE FILE FOR GATK
#samtools faidx polished_assembly.fasta


##########################################################################
########################### SET USEFUL PATHS #############################
##########################################################################

ref=/proj/julianog/refs/AbHGAP_George/polished_assembly.fasta
reads=/proj/julianog/sequence_reads/htsf_seq_backups/2014_11_24_Hajime_acineto_illumina
picard=/nas02/apps/picard-1.88/picard-tools-1.88
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
#dir=aln/ # This sets the working directory
dir=aln/ # This sets the working directory

##########################################################################
###################### SAMPLE ALIGNMENT & CLEANING #######################
##########################################################################

#for file in `cat names/filenames.txt`
#do

### ALIGN PAIRED-END READS WITH BWA_MEM
#bwa mem -M \
#	-t 2 \
#	-v 2 \
#	-R "@RG\tID:$file\tPL:illumina\tLB:$file\tSM:$file" \
#	$ref \
#	$reads/$file\_R1_001.fastq.gz \
#	$reads/$file\_R2_001.fastq.gz \
#	> $dir/$file.sam
##		 -M marks shorter split hits as secondary
##		 -t indicates number of threads
##		 -v 2 is verbosity ... warnings and errors only

##### POTENTIALLY MORE SENSITIVE ALIGNMENT
##bwa mem -M \
##	-t 2 \
##	-v 2 \
##	-A 3 \
##	-L 15 \
##	-U 9 \
##	-T 75 \
##	-k 19 \
##	-w 100 \
##	-d 100 \
##	-r 1.5 \
##	-c 10000 \
##	-B 4 \
##	-O 6 \
##	-E 1 \
##	-R "@RG\tID:$file\tPL:illumina\tLB:$file\tSM:$file" \
##	$ref \
##	$reads/$file\_R1_001.fastq.gz \
##	$reads/$file\_R2_001.fastq.gz \
##	> $dir/$file.sam
##		# -M marks shorter split hits as secondary
##		# -t indicates number of threads
##		# -v 2 is verbosity ... warnings and errors only

### SORT SAM FILE AND OUTPUT AS BAM
#java -jar $picard/SortSam.jar \
#	I=$dir/$file.sam \
#	O=$dir/$file.sorted.bam \
#	SO=coordinate

### MARK DUPLICATES
#java -jar $picard/MarkDuplicates.jar \
#	I=$dir/$file.sorted.bam \
#	O=$dir/$file.dedup.bam \
#	METRICS_FILE=$dir/$file.dedup.metrics \
#	REMOVE_DUPLICATES=True

### INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar $picard/BuildBamIndex.jar \
#	INPUT=$dir/$file.dedup.bam

### IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar $gatk \
#	-T RealignerTargetCreator \
#	-R $ref \
#	-I $dir/$file.dedup.bam \
#	-o $dir/$file.realigner.intervals

### PERFORM THE ACTUAL REALIGNMENT
#java -jar $gatk \
#	-T IndelRealigner \
#	-R $ref \
#	-I $dir/$file.dedup.bam \
#	-targetIntervals $dir/$file.realigner.intervals \
#	-o $dir/$file.realn.bam

#rm aln/*.sam
#rm aln/*dedup*
#rm aln/*sorted*
#rm aln/*intervals

### CALCULATE COVERAGE
##bedtools genomecov -ibam $dir/$file.realn.bam -max 10 | grep genome > coverage/aln_vs_HGAP/$file.cov10
##tail -1 coverage/aln_vs_HGAP/$file.cov10 | cut -f 5 >> coverage/aln_vs_HGAP/genome10.txt

#done

##########################################################################
############################ VARIANT CALLING #############################
##########################################################################

### MULTIPLE-SAMPLE VARIANT CALLING (UG IS GATK'S CALLER FOR NON-DIPLOID)
#java -jar $gatk \
#	-T UnifiedGenotyper \
#	-R $ref \
#	-I names/good42bamnames.list \
#	-o variants/good42_UG.vcf \
#	-ploidy 1 \
#	-nt 1

### TRYING OUT HAPLOTYPE CALLER
#java -jar $gatk \
#	-T HaplotypeCaller \
#	-R $ref \
#	-I names/good42bamnames.list \
#	-o variants/good42_HC.vcf \
#	-ploidy 1
#		# gatk.intervals includes just the chromosomes and mitochondria
#		# HC does not support -nt

##########################################################################
########################### VARIANT FILTERING ############################
##########################################################################

## FILTER BY DEPTH IN PERCENTAGE OF SAMPLES
java -Xmx2g -jar $gatk \
	-T CoveredByNSamplesSites \
	-R $ref \
	-V variants/good42_HC.vcf \
	-out variants/5xAT100%.intervals \
	-minCov 5 \
	-percentage 0.9999999
		# Output interval file contains sites that passed

## FILTER VCF BY DEPTH ACROSS SAMPLES
java -jar $gatk \
	-T VariantFiltration \
	-R $ref \
	-V variants/good42_HC.vcf \
	-L variants/5xAT100%.intervals \
	--logging_level ERROR \
	-o variants/good42_HC_5xAT100%.vcf
		# --logging_level ERROR suppresses any unwanted messages

## NOW NEED TO RUN vcfFilter.Rmd
## AND USE THOSE PLOTS TO DETERMINE QUAL VALUES TO CUT AT 

### FILTER VCF BY NEAFSEY PARALOGS, TANDEM REPEATS, SUBTELOMERES, AND QUALITY SCORES
#java -jar $gatk \
#	-T VariantFiltration \
#	-R $ref \
#	-V variants/good42.5xAT100%.vcf \
#	--filterExpression "QD < 20.0" \
#	--filterName "QD" \
#	--filterExpression "MQ < 55.0" \
#	--filterName "MQ" \
#	--filterExpression "FS > 5.0" \
#	--filterName "FS" \
#	--filterExpression "MQRankSum < -5.0" \
#	--filterName "MQRankSum" \
#	--filterExpression "ReadPosRankSum < -5.0" \
#	--filterName "ReadPosRankSum" \
#	--logging_level ERROR \
#	-o variants/good42.qual.vcf
#		# --logging_level ERROR suppresses any unwanted messages
#		# The three .intervals files contain intervals that should be excluded

##########################################################################
####################### FORMATTING FOR OUTBREAKER ########################
##########################################################################

#for name in `cat names/good42.txt`
#do

### SPLIT VCF INTO INDIVIDUAL VCFs
#java -Xmx2g -jar $gatk \
#	-T SelectVariants \
#	-R $ref \
#	--variant variants/good42.qual.vcf \
#	-sn $name \
#	-o variants/indivs/$name.vcf

### REMOVE ALL NON-ENTRIES IN INDIVIDUAL FILES (GT=0)
#grep -vP "PL\t0" variants/indivs/$name.vcf > variants/indivs/$name.sans0.vcf

### VCF TO FASTA
#java -Xmx2g -jar $gatk \
#	-T FastaAlternateReferenceMaker \
#	-R $ref \
#	--variant variants/indivs/$name.sans0.vcf \
#	-o variants/indivs/$name.fa
#		## --rawOnelineSeq prints only sequence

### CREATING THE MULTIFASTA
#echo ">"$name >> variants/good42.fa
#grep -v ">" variants/indivs/$name.fa >> variants/good42.fa

#done


##########################################################################
############################## EXTRA TOOLS ###############################
##########################################################################

#for bam in `cat names/needDeeperCov.txt`
#do

### MAKE IGV BATCH FILES
#echo "## Batch script to take pictures of the glta sequence
#new
#genome HGAP_2Cells_5Feb2015
#snapshotDirectory /proj/julianog/users/ChristianP/acinetoWGS/mlst-cleanup-A3
#load /proj/julianog/users/ChristianP/acinetoWGS/$dir/$bam.realn.bam
#goto unitig_31|quiver:232,505-233,139
#snapshot $bam.fusa.png
#goto unitig_30|quiver:1,996,986-1,997,470
#snapshot $bam.glta.png
#goto unitig_30|quiver:932,307-932,605
#snapshot $bam.pyrg.png
#goto unitig_30|quiver:1,007,379-1,007,752
#snapshot $bam.reca.png
#goto unitig_29|quiver:54,233-54,564
#snapshot $bam.rplb.png
#goto unitig_28|quiver:27,764-28,221
#snapshot $bam.rpob.png
#exit" > mlst/$bam.batch

### RUN THE BATCH SCRIPT IN IGV
#igv -b mlst/$bam.batch

#done


### GATK DEPTH OF COVERAGE CALCUALTOR
#java -Xmx10g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
#	-T DepthOfCoverage \
#	-R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta \
#	-I bamnames.list \
#	-o coverage/allExons.cov \
#	-geneList coverage/PvSal1_v10.0_exons.refseq \
#	-L coverage/PvSal1_v10.0_exons.intervals \
#	-omitBaseOutput \
#	--minMappingQuality 20 \
#	--minBaseQuality 20
#		#-omitBaseOutput gets rid of the large by-sample-by-base output file
#		# Apparently, we can provide a refseq file of features in the genome for site-by-site analysis
#		# http://gatkforums.broadinstitute.org/discussion/1329/using-refseq-data

## COUNT READS
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T CountReads -R $ref -I $dir/$name.merged.bam -rf MappingQualityZero

## SORT SAM FILE AND OUTPUT AS BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=$dir/$name-lane1.sam O=$dir/$name-lane1.sorted.bam SO=coordinate
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=$dir/$name-lane2.sam O=$dir/$name-lane2.sorted.bam SO=coordinate

## INDEX BAM FILE
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=$dir/1737Pv.sorted.bam

## VALIDATE VCF FORMAT FOR GATK
#java -Xmx2g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $ref -T ValidateVariants --validationTypeToExclude ALL --variant plasmoDB_vivax_snps.vcf

## COMPARE VCF FILES
#vcftools --vcf bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bt2.vcf --diff bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bwa.vcf --out bwa_vs_bt2/compare.txt
