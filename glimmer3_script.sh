#!/bin/bash

# /usr/bagel3/glimmer3_script.sh ./my_session/query.fna ./my_session/glimmer /usr/local/glimmer3.02/bin


genome=$1
tag=$2
glimmerpath=$3

echo genome = $genome
echo $glimmerpath

$glimmerpath/long-orfs -n -t 1.15 $genome $tag.longorfs > /dev/null 
# Extract the training sequences from the genome file
$glimmerpath/extract -t $genome $tag.longorfs > $tag.train 
# Build the icm from the training sequences
$glimmerpath/build-icm -r $tag.icm < $tag.train 
# Run Glimmer
$glimmerpath/glimmer3 -o50 -g90 -t30 --linear $genome $tag.icm $tag > /dev/null 

