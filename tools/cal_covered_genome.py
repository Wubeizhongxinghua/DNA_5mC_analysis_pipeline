import argparse
import sys
import subprocess
import re
import time
import os
import pandas as pd


def filtered_reads(file1,outname_prx):
	t1 = time.time()
	output_1 = outname_prx + "_covered.info"
	filex = open(file1)
	sum_0 = 0
	sum_N = 0
	sum_covered = 0
	for line in filex:
		line2  = line.strip().split("\t")
		cover_x = int(line2[2])
		if cover_x !=0:
			sum_N+=1
			sum_covered += cover_x
		else:
			sum_0 += 1
	xx1 = round(sum_N/(sum_N+sum_0)*100,2)
	print(sum_covered,sum_N+sum_0)
	aver = round(sum_covered/(sum_N),2)
	pd_data = pd.DataFrame({'Covered_sites':[sum_N],'Non-Covered_sites':[sum_0],'Covered_genome':[xx1],'Average_depth':[aver]})
	pd_data.to_csv(output_1,sep="\t",index=False)
	print(pd_data)






if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="get mutation rates for m6A_new")
	parser.add_argument("-file1", "--file1",nargs="?", type=str, default=sys.stdin, help="file1")
	parser.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")

	options = parser.parse_args()

	subprocess.call("samtools depth -a " + options.file1+" > " + options.outname_prx+".depth", shell=True)
	print("************************begin calculating*****************************")
	filtered_reads(options.outname_prx+".depth",options.outname_prx)
	os.system(f"rm -rf {options.outname_prx}.depth")
