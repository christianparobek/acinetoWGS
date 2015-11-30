## Count the nubmer of diffs between two VCF input files
## To respond to Reviewer #2's comment from AAC submission
## Started 25 November 2015
## This is the template to count differences
## Runs the diffCounter.sh engine



bash diffCounter.sh 03 11
bash diffCounter.sh 03 16
bash diffCounter.sh 03 14
bash diffCounter.sh 03 08
bash diffCounter.sh 16 08
bash diffCounter.sh 14 07
bash diffCounter.sh 14 10
bash diffCounter.sh 08 05
bash diffCounter.sh 08 01
bash diffCounter.sh 08 18
bash diffCounter.sh 08 09
bash diffCounter.sh 01 02
bash diffCounter.sh 18 15
bash diffCounter.sh 09 18
bash diffCounter.sh 09 15
bash diffCounter.sh 09 17
bash diffCounter.sh 08 17
bash diffCounter.sh 15 25
bash diffCounter.sh 15 36
bash diffCounter.sh 20 19
bash diffCounter.sh 20 22
bash diffCounter.sh 20 26
bash diffCounter.sh 20 21
bash diffCounter.sh 20 27
bash diffCounter.sh 20 23
bash diffCounter.sh 19 22
bash diffCounter.sh 19 21
bash diffCounter.sh 19 26
bash diffCounter.sh 19 23
bash diffCounter.sh 19 27
bash diffCounter.sh 22 24
bash diffCounter.sh 22 30
bash diffCounter.sh 22 28
bash diffCounter.sh 28 33
bash diffCounter.sh 33 29
bash diffCounter.sh 33 31
bash diffCounter.sh 30 37
bash diffCounter.sh 30 35
bash diffCounter.sh 35 43
bash diffCounter.sh 43 45
bash diffCounter.sh 43 46
bash diffCounter.sh 43 41
bash diffCounter.sh 46 48
bash diffCounter.sh 46 39
bash diffCounter.sh 41 40
bash diffCounter.sh 41 38
bash diffCounter.sh 41 34
bash diffCounter.sh 38 40
bash diffCounter.sh 38 34


for i in {01..48}
do

bash diffCounter.sh 32 $i

done

