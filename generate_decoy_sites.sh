#!/bin/bash

# This is the code for generating decoy RNA A-to-I editing sites
# 
# Required files and modules:
# original_REDI_sites with format "chr1:11111111"
# corresponding reference genome
# samtools

module load samtools

coor_plus_one() {

        chr=$(echo $1|cut -d":" -f1)
        coor=$(echo $1|cut -d":" -f2)

        coor=$((coor+1))

        new="$chr:$coor"

        echo $new

}

generate_pattern() {

	chr=$(echo $1|cut -d":" -f1)
        coor=$(echo $1|cut -d":" -f2)

        new="$chr:$coor-$coor"

	echo $new
}



generate_new() {

	FILE1="original_REDI_sites"
	FILE2="decoy_sites"
	STRING="$1"
	BASE="C"
	a=$(grep -x "$STRING" "$FILE1")
	b=$(grep -x "$STRING" "$FILE2")

        while [ ! -z "$a" ] || [ ! -z "$b" ] || [[ "$BASE" =~ ^(c|g|C|G)$ ]]
	do
                STRING=$(coor_plus_one $STRING)
		PATTERN=$(generate_pattern $STRING)
		BASE=$(samtools faidx hg19.fasta $PATTERN | tail -1)
		a=$(grep -x "$STRING" "$FILE1")
		b=$(grep -x "$STRING" "$FILE2")
        done
	echo $STRING
	echo $BASE
}

touch decoy_sites

for i in $(cat original_REDI_sites)
do
	generate_new $i >> decoy_sites
done

grep chr decoy_sites > coor
grep -v chr decoy_sites > base
sed 's/\(A\|a\)/+/g; s/\(T\|t\)/-/g' base > strand

awk -v FS=':' -v OFS='\t' '{print $1,$2}' coor > tmp1
paste tmp1 strand > decoy_sites

rm coor base strand tmp1




