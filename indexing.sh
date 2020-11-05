#!/bin/bash

for i in *bam_sorter.bam

do

echo "Indexing: "$i        

samtools index $i $i".bai"

done
