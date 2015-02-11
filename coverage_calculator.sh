### Coverage Calculation pipeline for Acinetobacter baumannii whole-genome sequencing
### Started December 19, 2014

for name in `cat bamnames.list | head -2`
do

echo $name - 
bedtools genomecov -ibam $name -max 10 | grep "genome	10"

#| grep "genome 10" >> coverage/genome10.txt

done


genome	10
