#!/bin/bash
for spl in HepG2-WGA-GLORI-R2_L1 HepG2-WGA-untreated_L4
do
	pkurun-cnlong 1 8 python3 cal_covered_genome.py -file1 ${spl}/${spl}_merged.rmdup.sorted.bam -outname_prx ${spl}
	sleep 1
	pkurun-cnlong 1 8 python3 cal_covered_genome.py -file1 ${spl}/${spl}_hg38_merged.rmdup.sorted.bam -outname_prx ${spl}_hg38
	sleep 1
done
