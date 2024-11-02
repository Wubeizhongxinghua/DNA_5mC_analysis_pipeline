rule trim_galore:
	input:
		fastq1 = "data/{dataset}/{sample_id}_1.fq.gz",
		fastq2 = "data/{dataset}/{sample_id}_2.fq.gz"
	output:
		fastq1 = temp("output/{dataset}/trim/{sample_id}_1_val_1.fq.gz"),
		fastq2 = temp("output/{dataset}/trim/{sample_id}_2_val_2.fq.gz")
	threads: 8
	shell:
		"""
		trim_galore -j 20 -q 30 --phred33 --length 35 --stringency 3 --paired -o output/{wildcards.dataset}/trim/ {input.fastq1} {input.fastq2}
		"""
# todo: add multiqc function?
rule fastqc:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		report1 = "output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html",
		report2 = "output/{dataset}/trim/{sample_id}_2_val_2_fastqc.html"
	threads: 20
	shell:
		"""
		fastqc -o output/{wildcards.dataset}/trim/ -t {threads} {input.fastq1} {input.fastq2}
		"""
# alignment with bwameth, you need to check whether your toolshed package has been installed.
rule bwameth:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		bam = temp("output/{dataset}/bwameth_results/{sample_id}.bam")
	log:
		"output/{dataset}/log/bwameth_results/{sample_id}.log"
	threads: 20
	shell:
		"""
		bwameth.py --reference {genome_dir}/{genome_fasta_file_name} -t {threads} {input.fastq1} {input.fastq2} 2> {log} | samtools view -h -b -o {output.bam}
		"""

rule sort:
	input:
		bam = "output/{dataset}/bwameth_results/{sample_id}.bam"
	output:
		bam = temp("output/{dataset}/bwameth_results/{sample_id}_sort.bam")
	log:
		"output/{dataset}/log/bwameth_results/{sample_id}_sort.log"
	threads: 20
	shell:
		"""
		samtools sort -@ 20 {input.bam} -o {output.bam}
		samtools index -@ 20 {output.bam}
		"""

# rmduplicate using picard.jar
rule rmduplicate:
	input:
		bam = "output/{dataset}/bwameth_results/{sample_id}_sort.bam"
	output:
		bam = temp("output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.deduplicated.bam")
	log:
		"output/{dataset}/log/bwameth_results/bam_rmdup/{sample_id}.log"
	threads: 20
	shell:
		"""
		java -jar tools/picard.jar MarkDuplicates -I {input.bam} -O {output.bam} -REMOVE_DUPLICATES true -METRICS_FILE output/{wildcards.dataset}/bwameth_results/bam_rmdup/{wildcards.sample_id}_rmdup.txt -CREATE_INDEX true -USE_JDK_DEFLATER true -USE_JDK_INFLATER true --TMP_DIR tmpdir
		"""

rule filter:
	input:
		bam = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.deduplicated.bam"
	output:
		bam = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam",
		bamstat = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam.stat"
	log:
		"output/{dataset}/log/bwameth_results/bam_rmdup/{sample_id}_clean.log"
	threads: 20
	shell:
		"""
		samtools view -h -@ 20 -q 60 -F 780 -f 2 {input.bam} | samtools sort -@ 20 -n | python3 tools/remove_3CH_readpairs.py - | samtools sort -@ 20 > {output.bam}
		samtools index -@ 20 {output.bam}
		samtools stat -@ 20 {output.bam} > {output.bamstat}
		"""

rule insert_size:
	input:
		bam = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam"
	output:
		metrics = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}_size_matrics.txt"
	threads: 3
	shell:
		"""
		java -jar tools/picard.jar CollectInsertSizeMetrics I={input.bam} O={output.metrics} H=output/{wildcards.dataset}/bwameth_results/bam_rmdup/{wildcards.sample_id}_size_matrics.pdf M=0.5
		"""

rule stat_readnum:
	input:
		bamstat = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam.stat"
	output:
		bamnum = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.num.txt"
	threads: 1
	shell:
		"""
		echo -e "{wildcards.sample_id}\t$(grep -e "SN\ssequences" {input.bamstat} | cut -d ':' -f 2 | grep -o -e "[0-9]*")" > {output.bamnum}
		"""

rule index:
	input:
		bamclean = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam"
	output:
		cleanbam_idx = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam.bai"
	threads: 20
	shell:
		"""
		samtools index -@ {threads} {input.bamclean}
		"""

rule cal_coverage:
	input:
		bam = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam"
	output:
		covinfo = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}_cover_genome_covered.info"
	threads: 5
	log:
		"output/{dataset}/log/bwameth_results/bam_rmdup/{sample_id}_cov.log"
	shell:
		"""
		python3 tools/cal_covered_genome.py -file1 {input.bam} -outname_prx output/{wildcards.dataset}/bwameth_results/bam_rmdup/{wildcards.sample_id}_cover_genome
		rm -rf output/{wildcards.dataset}/bwameth_results/bam_rmdup/{wildcards.sample_id}_cover_genome.depth
		"""


rule reads:
	input:
		bam = "output/{dataset}/bwameth_results/{sample_id}_sort.bam"
	output:
		bamread = "output/{dataset}/bwameth_results/{sample_id}_sort_reads.txt"
	threads: 20
	shell:
		"""
		samtools view -F 260 -@ {threads} -c {input.bam} > {output.bamread}
		"""


rule methyldackel:
	input:
		bam_idx = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam.bai",
		bam_rmdup = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}.clean.bam"
	output:
		bedgraph = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}_CpG.bedGraph"
	log:
		"output/{dataset}/log/bwameth_results/methyldackel/{sample_id}.log"
	threads: 20
	shell:
		"""
		bash tools/methyldackel.sh output/{wildcards.dataset}/bwameth_results/bam_rmdup/ {wildcards.sample_id} {genome_dir}/{genome_fasta_file_name} {threads}
		echo "Rerun"
		"""

rule total_beta:
	input:
		bedgraph = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}_CpG.bedGraph"
	output:
		total_beta = "output/{dataset}/bwameth_results/bam_rmdup/{sample_id}_totalbeta.txt"
	threads: 5
	run:
		import polars as pl
		bed = pl.read_csv(input.bedgraph, separator='\t', has_header=False, new_columns=['chr','start','end','beta','m','um'], comment_prefix='t')
		total_betavalue = bed.with_columns(
			(
				pl.col('m').sum()
				/
				(
					pl.col('m').sum()
					+
					pl.col('um').sum()
				)
			).alias('total_beta')
		)['total_beta'][0]
		with open(output.total_beta, 'w') as f:
			f.write(f'{total_betavalue}\n')


##################################### align to qc ##########################33
#
rule bismark_lamda:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		bam = "output/{dataset}/bismark_results_lamda/{sample_id}_1_val_1_bismark_bt2_PE_report.txt"
	log:
		"output/{dataset}/log/bismark_results_lamda/{sample_id}.log"
	threads: 20
	shell:
		"""
		bismark --bowtie2 -p {threads} -N 1 -L 30 {lamda_dir} -q -1 {input.fastq1} -2 {input.fastq2} -o output/{wildcards.dataset}/bismark_results_lamda/ > {log} 2>&1
		"""



rule bismark_pU:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		bam = "output/{dataset}/bismark_results_pU19/{sample_id}_1_val_1_bismark_bt2_PE_report.txt"
	log:
		"output/{dataset}/log/bismark_results_pU19/{sample_id}.log"
	threads: 20
	shell:
		"""
		bismark --bowtie2 -p {threads} -N 1 -L 30 {pU19_dir} -q -1 {input.fastq1} -2 {input.fastq2} -o output/{wildcards.dataset}/bismark_results_pU19/ > {log} 2>&1
		"""
