#!/bin/bash
#bams=("HepG2-GLORI-R1_merge_hg38_merged.sorted.bam" "HepG2-GLORI-R1_merge_hg38_merged.rmdup.sorted.bam" "HepG2-GLORI-R1_merge_hg38_merged_outpeak.rmdup.sorted.bam")
for sml in HepG2-WGA-GLORI-R2_L1 HepG2-WGA-untreated_L4
do
	for bam in `dir -1 ./${sml}/${sml}*merged.sorted.bam`
	do
		pkurun-cnlong 1 5 "samtools view -@ 5 -c ${bam} > ${bam}_reads.txt"
		sleep 1
	done
done
