'''
Copyright (c) 2024-03-26 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs3/chengqiyi_pkuhpc/limingyang/model/tools/check_read_paired.py

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
'''
from rich import print, pretty
from rich.traceback import install
pretty.install()
install(show_locals=True)
import pysam
import sys

def read_pair_generator(bam):
	"""
	Generate read pairs from a BAM file. Pairs are identified based on their query name (QNAME).
	Assumes BAM is sorted by name.
	"""
	# Dictionary to temporarily store reads until their pair is found
	read_pairs = {}

	for read in bam:
		if read.query_name in read_pairs:
			yield read, read_pairs.pop(read.query_name)
		else:
			read_pairs[read.query_name] = read

	# Yield remaining reads if their pair wasn't found
	for read in read_pairs.values():
		yield read, None

def contains_three_or_more_special_sequences(read, flag):
	"""Check if a read contains 3 or more instances of CC, CT, or CA."""
	total_count = 10 #initialize
	if flag in [147, 99]: # R1 +
		sequences = ['CC', 'CT', 'CA']
		total_count = sum(read.query_sequence.count(seq) for seq in sequences)
	elif flag in [83, 163]: # R1 -
		sequences_rev = ['GG', 'AG', 'TG']
		total_count = sum(read.query_sequence.count(seq) for seq in sequences_rev)
	return total_count >= 3



# Using stdin and stdout for input and output
with pysam.AlignmentFile("-", "rb", check_sq=False) as infile, pysam.AlignmentFile("-", "wb", template=infile) as outfile:
	for read1, read2 in read_pair_generator(infile): # gen paried reads
		if read1 and read2:
			if not contains_three_or_more_special_sequences(read1, read1.flag) and not contains_three_or_more_special_sequences(read2, read2.flag): #read1 and 2 passed the filter
				outfile.write(read1)
				outfile.write(read2)

# Important: Flush the output to ensure all data is written out
sys.stdout.flush()

