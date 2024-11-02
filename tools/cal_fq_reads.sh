#!/bin/bash
for i in `dir -1 *_1.fq.gz`
do
	pkurun-cns 1 4 "zcat ${i} | grep -c "^+$" > ${i%%_1.fq.gz}_readnum.txt"
	sleep 1
done
