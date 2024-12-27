'''
Copyright (c) 2024-03-26 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs3/chengqiyi_pkuhpc/limingyang/model/tools/remove_3CH_reads.py

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

def contains_three_or_more_special_sequences(read, is_reverse):
	"""Check if a read contains 3 or more instances of CC, CT, or CA."""
	if not is_reverse: # +
		sequences = ['CC', 'CT', 'CA']
		total_count = sum(read.query_sequence.count(seq) for seq in sequences)
	else: # -
		sequences_rev = ['GG', 'GA', 'GT']
		total_count = sum(read.query_sequence.count(seq) for seq in sequences_rev) 
	return total_count >= 3

# 使用 stdin 和 stdout 作为输入和输出
with pysam.AlignmentFile("-", "rb", check_sq=False, check_header=False) as infile, \
	pysam.AlignmentFile("-", "wb", template=infile) as outfile:
	for read in infile:
		if not contains_three_or_more_special_sequences(read, read.is_reverse):
			outfile.write(read)

