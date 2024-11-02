## **DNA 5mC Methylation Analysis Pipeline**

This repo. provides a fastq-to-final-result straightforward pipeline based on `snakemake` to conduct convenient 5mC analysis suitable for **paired-end** BS-seq data or other libraries with similar principle (C to T, 5mC to C).

## **Deployment and Setup**

1. Clone this repo.: `git clone https://github.com/Wubeizhongxinghua/DNA_5mC_analysis_pipeline.git`

2. Navigate to the repo. dir.: `cd DNA_5mC_analysis_pipeline`

3. Put your dataset directory to `./data`, where the name of dataset directory is your dataset name. Inside the directory, all data files are named like `samplename_1.fq.gz` and `samplename_2.fq.gz`

```
.
├── data
│   └── example_dataset
│       ├── example_sample_1.fq.gz
│       └── example_sample_2.fq.gz
```

Here, the dataset name is `example_dataset`, where there is one sample called `example_sample`.

4. Collect all the sample names into a .txt file, named using the dataset name. Put the file into `./sample_id`

```
├── sample_id
│   └── example_dataset.txt
```

```shell
cat sample_id/example_dataset.txt

# return: example_sample
```

The format of the .txt file should satisfy:
- Each row is a sample name.
- You can add other metadata of the sample, like its group or sth. in the sample row with **delim \t**. Like:
```
example_sample	cancer	age35
example_sample2	healthy	age32
```
- There should not be any header of the file.

5. Edit the `methylation_BS_EMseq.snake` for settings:

- Edit the `genome_dir` to the dir. of your genome.
- Edit the `genome_fasta_file_name` to the genome fasta file name, which must be inside the `genome_dir` dir.
- Edit the `lamda_dir` to the lambda DNA dir. name, where you should set the `bismark` index.
- Edit the `pUC19_dir` to the pUC19 DNA dir. name, where you should set the `bismark` index.
- Edit the list `datasets_all`, `datasets_nopUC19`, `datasets_nolambda`, `datasets_noBSQC`, where you shuold put the dataset name inside.
	- If you want the all analysis pipeline conducted for the dataset, add the dataset name into `datasets_all`.
	- If you want the analysis pipeline conducted for the dataset **except alignment to pUC19**, add the dataset name into `datasets_nopUC19`.
	- If you want the all analysis pipeline conducted for the dataset **except alignment to lambda DNA**, add the name into `datasets_nolambda`.
	- If you want the all analysis pipeline conducted for the dataset **except alignment to lambda DNA and pUC19**, add the name into `datasets_noBSQC`.

6. Install all python requirements and softwares:

- Softwares
	- `snakmeake` (Only tested using snakemake 7.32.3. Not sure whether snakemake 8 is suitable or not.)
	- `trim_galore`
	- `fastqc`
	- `bwameth`
	- `bismark`
	- `samtools`
	- `methylDackel`
	
- Python modules
	- `polars`
	- `pandas`
	- `pysam`
	- `rich`
	- `toolshed` (required by [bwameth](https://github.com/brentp/bwa-meth))

7. Test your pipe line using:

```shell
snakemake -s methylation_BS_EMseq.snake --dry-run --rerun-triggers mtime -pr 
```

8. If everything works, then you can run the pipeline.
