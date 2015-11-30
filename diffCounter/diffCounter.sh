## Count the nubmer of diffs between two VCF input files
## To respond to Reviewer #2's comment from AAC submission
## Started 25 November 2015
## This is the engine to count the number of differences
## Nests inside a diffCounterTemplate.sh script

vcfdir=/proj/julianog/users/ChristianP/acinetoWGS/variants/indivs_UG/
#ref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
#gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
#vcfdir=/proj/julianog/users/ChristianP/amaAbs/variants/split_indivs_HC

###########################################
##### GET VARIANTS THAT AREN'T SHRAED #####
###########################################


bedtools intersect -v -a $vcfdir/CLIN_A$1*vcf -b $vcfdir/CLIN_A$2*vcf > in_A_not_in_B.vcf
		# -v reports only entries in A w no overlaps in B

bedtools intersect -v -a $vcfdir/CLIN_A$2*vcf -b $vcfdir/CLIN_A$1*vcf > in_B_not_in_A.vcf
		# -v reports only entries in A w no overlaps in B

echo $1 $2
cat in_A_not_in_B.vcf in_B_not_in_A.vcf | wc -l
