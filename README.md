# RNAseq-pipeline (Paired-end)
Author: Jeanne-Marie McPherson

# Quick Start:
## Paired-end RNA-seq
### Clone pipeline (Current stable version)
git clone https://github.com/mckaylabunc/RNAseq-pipeline.git 

Make sure all necessary files are present: 

        --sampleSheet.csv
        --runSnakemake.sh
        --snakefile_RNA_snakemake.py
        --config_RNA.json
        --slurmConfig

### Create sampleSheet_EXAMPLE_RNAseq.csv with descriptive columns of data.
_Note: You need separate rows for fastq_R1 and fastq_R2 for paired-end sequencing data._

_You can include further descriptive columns that may be useful for downstream DESeq2 analysis._

	htsfFile	 		fastq			sampleName		
	path/to/mySample1_R1.fastq.gz	mySample1_R1.fastq.gz	mySample1_rep#	
	path/to/mySample1_R2.fastq.gz	mySample1_R2.fastq.gz 	mySample1_rep#
	path/to/mySample2_R1.fastq.gz	mySample2_R1.fastq.gz	mySample2_rep#	
	path/to/mySample2_R2.fastq.gz	mySample2_R2.fastq.gz 	mySample2_rep#

### Set parameters in config_RNA.json
1. References correct sampleSheet.

2. Edit readLen as needed.

### Edit slurmConfig.json to configure default parameters if necessary.

### Edit parameters in `rule STARindex` based on readLen & genome you want to align to.
The flag `--sjdbOverhang` will change based on your read length. Make sure to set to readLen-1. For example, if your reads are 75bp, you should have `--sjdbOverhang 74`.

If you plan to align reads to genome other than dm6 you will need to change the input files in `rule STARindex`.

### Deploy submission with runSnakemake.sh
Command to run pipeline within the RNAseq-pipeline directory: 

        $ bash runSnakemake.sh snakefile_RNA_snakemake.py 
        
Command to run pipeline in the background:

        $ sbatch --mem=256MB -e snakeErr.%J -o snakeOut.%J --wrap='bash runSnakemake.sh snakefile_RNA_snakemake.py'


### Other factors to keep in mind:
Is your RNAseq library prep kit stranded or reverse stranded? This will change flag in `rule BamtoBigWig`. `--filterRNAstrand reverse` is assuming a reverse stranded library prep kit. Use `--filterRNAstrand forward` if using a stranded library prep kit.

# Description of Pipeline Steps
This pipeline in implemented in Snakemake, a workflow manager for handling job dependencies and submissions to HPC clusters. Snakemake will automate the job submission process for each sample individually and will therefore often run some of these steps seemingly out of order. Below is described the general flow of information through the pipeline. In some instances many jobs are handled by snakemake to accomplish a single defined "step" below.

The rules found in the Snakefile for this pipeline are written generally in the order in which they happen from top to bottom. See there for more detail.

1. Combine technical replicates.
	- Autodetects technical replicates and pools reads for alignment.
2. Generate STAR index.
3. Adapter trimming using bbduk.
	- soft-clips illumina adapters from reads.
4. Alignment with STAR.
5. Sort Bam file & convert to bed format with bedtools.
6. Make coverage files (bigwig format).
	- Generates files for combined, forward, and reverse strands. Used for visualization in IGV.
7. Read summarization with featureCounts.
	- Generates .txt file that can be used as input for DESeq2 analysis. 
8. Quality control.
	- Fastqc, multiqc.
	- Library complexity with Picard to determine % duplicates.

# Troubleshooting
When developing for a new compute environment, running `snakemake -n -p` will run in "dry-run" mode, which can be useful for debugging issues with configuration before running the data processing steps.

Consider setting `--rerun-incomplete` in the Snakemake call in runSnakemake.sh when initially implementing the pipeline on a new compute environment or testing new parameters.

### Common Error Messages 
If you get `Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory...`, then run this command to unlock the directory: 
	
	$ snakemake --snakefile snakefile_RNA_snakemake.py --unlock --cores 3 --rerun-incomplete
	
If you get `error: *** JOB 28546694 ON c0933 CANCELLED AT ___ DUE TO TIME LIMIT ***` then extend time limits in `slurmConfig` file.
