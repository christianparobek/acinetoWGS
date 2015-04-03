### Started April 1 2015
### To calculate differences between VCF files

good42names=/proj/julianog/users/ChristianP/acinetoWGS/names/good42.txt
vcfdir=/proj/julianog/users/ChristianP/acinetoWGS/variants/indivs_UG/

##############################
#### CALCULATE DIFFERENCES ###
##############################
for name1 in `cat $good42names`
do

# subset the name to something reasonable
substr1=`echo $name1 | cut -d '_' -f 2`

	for name2 in `cat $good42names`
	do
	
	# subset the name to something reasonable
	substr2=`echo $name2 | cut -d '_' -f 2`
	
	# do the comparison
	vcftools --vcf $vcfdir$name1.sans0.vcf --diff $vcfdir$name2.sans0.vcf
	
	# move the files to a subdir
	mv out.diff.indv_in_files comparisons/$substr1-$substr2.indiv_in_files
	mv out.diff.sites_in_files comparisons/$substr1-$substr2.sites_in_files
	mv out.log comparisons/$substr1-$substr2.log
	done

done



#############################
### ADD DIFFS TO PFGE RES ###
#############################

echo -e "Comparison\tPFGE_status\tNb_variants" > variants_vs_PFGE_results.txt
	# make a blank file for the results
	# the -e lets us echo a tab

for line in `cat PFGE_results.txt`
do

comparo=`echo $line | cut -d '_' -f 1`
	# get the comparison identifier
PFGE=`echo $line | cut -d '_' -f 2`
	# get the PFGE status
vcf1diffs=`cat comparisons/$comparo.log | grep "SNPs only in main file" | cut -d " " -f 2`
	# calculate number of variants unique to the first vcf file
vcf2diffs=`cat comparisons/$comparo.log | grep "SNPs only in second file" | cut -d " " -f 2`
	# calculate number of variants unique to the second vcf file
totaldiffs=`expr $vcf1diffs + $vcf2diffs`
	# sum the total number of diffs

echo -e "$comparo\t$PFGE\t$totaldiffs" >> variants_vs_PFGE_results.txt
	# print results to the file

done
