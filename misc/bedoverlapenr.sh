#!/bin/bash

file1=$1; # the peaks file
file2=$2; # the "context" file (promoters, genes, etc)	
genomefile=$3; # genome file
bs=$4; # number of bootstraps
#echo $file1 $file2 $genomefile $bs;
#
size1=$(gawk '{val+=($3-$2)}END{print val}' < $file1); # Raw coverage of Peak file
gsize=$(gawk '{sum+=$2}END{print sum}' < $genomefile); # Size of Genome
ratio1=$(gawk -v gs=$gsize '{val+=($3-$2)}END{print val/gs}' < $file2); # Coverage of context file
# Coverage of intersection 
ratio=$(intersectBed -wao -a $file1 -b $file2 | gawk -v size=$size1 '{if ($NF!=0) {val+=$NF}}END{print (val/size)}');
#
# Key Calculation: Enrichment of Peak/Context Overlap vs Expected Ratio
enr=$(intersectBed -wao -a $file1 -b $file2 | gawk -v size=$size1 -v ratio=$ratio1 '{if ($NF!=0) {val+=$NF}}END{print (val/size)/ratio}');
#
echo "Enrichment of overlap=" $enr
# in case enrichment is =0 skip bootstraps by setting bs=0
if [ $enr ==  0 ] 
	then
	bs=0
fi
# else case ($enr>0)
echo "Calculating significance for "$bs" permutations. (Grab a coffee, this may take a while...)"
# clean up of temp file
rm -f temp.out;
# Initializing bootstraps
	for ((i=1;i<=$bs; i++)); 
	do
		# creating random permutation of Peaks file
		# Results stored in temp file
		shuffleBed -i $file1 -g $genomefile > test.bed; intersectBed -wao -a test.bed -b $file2 | gawk -v size=$size1 -v ratio=$ratio1 '{if ($NF!=0) {val+=$NF}}END{print (val/size)/ratio}' >> temp.out; 
		echo $i | gawk '{if ($1%100==0){print "permutations="$1 }}'; 
	done
# calculating bootstrap p-value (number of times en is > or < than bootstrap repeats)
pvalue=$(gawk -v en=$enr '{if (((en>1) && ($1>=en)) || ((en<1) && ($1<=en))) {val++}}END{print val/NR}' < temp.out | gawk -v lim=$bs '{if ($1==0) {print 1/lim} if ($1!=0) {print $1}}'); 
	# in case enrichment is =0 p.value=1
	if [ $enr ==  0 ] 
	then
	pvalue=1
fi
	echo "File1= "$file1 "File2= "$file2 "Ratio1= "$ratio1 "Ratio= "$ratio "Enrichment= "$enr "p-value= "$pvalue >> results
	echo "File1= "$file1 "File2= "$file2 "Ratio1= "$ratio1 "Ratio= "$ratio "Enrichment= "$enr "p-value= "$pvalue
	echo "finished run of" $file1 $file2
