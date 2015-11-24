### Coverage Calculation pipeline for Acinetobacter baumannii whole-genome sequencing
### Started December 19, 2014

for name in `cat names/bamnames.list`
do

echo $name - >> coverage/genome10.txt
bedtools genomecov -ibam $name -max 10 | grep "genome	10" >> coverage/genome10.txt
bedtools genomecov -ibam $name -max 5 | grep "genome	5" >> coverage/genome10.txt
#| grep "genome 10" >> coverage/genome10.txt

done



#for name in `cat names/bamnames.list`
#do echo $name -
#bedtools genomecov -ibam $name -max 10 | grep "genome\t5"
#done
