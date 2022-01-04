## Pipeline made by JEANNE-MARIE MCPHERSON.
## Edited September 2021.

import pandas as pd
import os
import shutil
import re
import glob
import sys
import getpass
import argparse

## --------------------------------------------------------------------------------------##
## Sample file config
## --------------------------------------------------------------------------------------##

## Configuration file	
configfile: 'config_RNA.json'

sampleDF = pd.read_csv(config['sampleSheetPath'], comment = '#')

def getsampleReplicates(wildcards):
	readNumRegex = '_R{}'.format(wildcards.Num)

	sampleFilter = sampleDF [ sampleDF.sampleName == wildcards.sample ]
	fastqList = list(sampleFilter.fastq)
	fastqList = [ fastq for fastq in fastqList 
			if re.search(readNumRegex, fastq) ]
	fastqList = [ 'Fastq/{}'.format(fastq) for fastq in fastqList ]

	return(fastqList)

sampleList =  set(sampleDF.sampleName)

## --------------------------------------------------------------------------------------##
## Snakemake rules
## --------------------------------------------------------------------------------------##

localrules: all

rule all:
	input:
		expand('Fastq/{fastq}', fastq = list(sampleDF.fastq)),
		expand('Fastq/{sample}_R{Num}.fastq.gz', sample = sampleList, Num = ['1', '2']),
		expand('Fastq/{sample}_R{Num}_trim.fastq.gz', sample = sampleList, Num = ['1', '2']),
		expand('Fastqc/{sample}_R{Num}_fastqc.html', sample = sampleList, Num = ['1', '2']),
		expand('Fastqc_trimmed/{sample}_R{Num}_trim_fastqc.html', sample = sampleList, Num = ['1', '2']),
		expand('libComplexity/{sample}_libraryComplexity.txt', sample = sampleList),
		expand('multiqc_data/multiqc_{statType}.txt', statType = ['fastqc', 'general_stats', 'sources']),
		expand('multiqc_report.html'),
		expand('Bed/{sample}_Aligned.sortedByCoord.out.bed', sample = sampleList),
		expand('Bam/{sample}_Aligned.sortedByCoord.out.bam', sample = sampleList),
		expand('BigWig/{sample}_Aligned.sortedByCoord.bw', sample = sampleList),
		expand('BigWig/{sample}_Aligned_fwd.sortedByCoord.bw', sample = sampleList),
		expand('BigWig/{sample}_Aligned_rev.sortedByCoord.bw', sample = sampleList),
		expand('featureCounts/counts_Aligned.txt'),
		expand('Star_Index/genomeParameters.txt')

rule moveFiles:
	input:
		lambda x: list(sampleDF.htsfFile)
	output:
		expand('Fastq/{fastq}', fastq = list(sampleDF.fastq))
	shell:
		"""
		inputList=({input})
		for ((i=0; i<${{#inputList[@]}}; i++)); 
		do
			cp ${{inputList[$i]}} Fastq/
		done
		"""

rule combinesampleReps:
	input:
		getsampleReplicates
	output:
		'Fastq/{sample}_R{Num,[12]}.fastq.gz'
	shell:
		"""
		cat {input} > {output}
		"""

## --------------------------------------------------------------------------------------##
## Build STAR Index
## --------------------------------------------------------------------------------------##
## NOTE: you will need to change --sjdbOverhang based on readLen, should be readLen-1 
rule STARindex:
       input:
               genome = config['index']['dm6']['fasta'],
               gtf = config['general']['annotation']
       output:
               'Star_Index/genomeParameters.txt'

       params:
               module = config['module']['starVer'],
               readLen = config['general']['readLen']
       shell:
               """
               module purge && module load {params.module}
               star --runThreadN 8 --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir Star_Index --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang 74
               """

## --------------------------------------------------------------------------------------##
## Adapter trimming
## --------------------------------------------------------------------------------------##
## Note: You may need to change the ktrim parameter based on the sequencer used. Illumina sequencing usually uses 3' adapters, so ktrim = r. If 5' adapter then ktrim = l. Also may need to change k= based on adapter length and minlen based on read length.
rule adapter_trim_reads:
	input:
		readOne = 'Fastq/{sample}_R1.fastq.gz',
		readTwo = 'Fastq/{sample}_R2.fastq.gz',
	output:
		readOne = 'Fastq/{sample}_R1_trim.fastq.gz',
		readTwo = 'Fastq/{sample}_R2_trim.fastq.gz',
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	params:
		module = config['module']['bbmapVer'],
		readLen = config['general']['readLen']
	shell:
		"""
		module purge && module load {params.module}
		bbduk.sh in1={input.readOne} in2={input.readTwo} out1={output.readOne} out2={output.readTwo} stats={output.adapterStats} ktrim=r k=23 minlen={params.readLen} ref=adapters rcomp=t tbo=t tpe=t hdist=1 mink=11 > {output.trimStats}
		"""	

## --------------------------------------------------------------------------------------##
## Quality Control 
## --------------------------------------------------------------------------------------##

## Fastqc, original reads
rule fastqc:
	input: 
		fastq = expand('Fastq/{sample}_R{Num}.fastq.gz', sample = sampleList, Num = ['1', '2'])
	output:
		'Fastqc/{sample}_R{Num}_fastqc.html'
	params:
		module = config['module']['fastqcVer']
	shell:
		"""
		module purge && module load {params.module}
		fastqc -o ./Fastqc/ -f fastq {input.fastq} 
		"""

## Fastqc, trimmed reads
rule fastqc_trimmed:
	input: 
		fastq = expand('Fastq/{sample}_R{Num}_trim.fastq.gz', sample = sampleList, Num = ['1', '2'])
	output:
		'Fastqc_trimmed/{sample}_R{Num}_trim_fastqc.html'
	params:
		module = config['module']['fastqcVer']
	shell:
		"""
		module purge && module load {params.module}
		fastqc -o ./Fastqc_trimmed/ -f fastq {input.fastq} 
		"""
## Multiqc
rule multiQC:
	input:
		qc = expand('Fastqc_trimmed/{sample}_R{Num}_trim_fastqc.html', sample = sampleList, Num = ['1', '2']), 
	output:
		fastqc = 'multiqc_data/multiqc_fastqc.txt',
		stats = 'multiqc_data/multiqc_general_stats.txt',
		sources = 'multiqc_data/multiqc_sources.txt',
		report = 'multiqc_report.html'
	params:
		module = config['module']['multiqcVer']
	shell:
		"""
		module purge && module load {params.module}
		multiqc -f . -o ./
		"""

## Determine library complexity (% duplicates)
rule libComplexity:
                input:
                        bam = 'Bam/{sample}_Aligned.sortedByCoord.out.bam'
                output:
                        libComplexity = 'libComplexity/{sample}_libraryComplexity.txt'
                params:
                        module = config['module']['picardVer'],
                        picardPath = config['module']['picardPath']
                shell:
                        """
                        module purge && module load {params.module}
                        java -jar {params.picardPath} EstimateLibraryComplexity INPUT= {input.bam} OUTPUT= {output.libComplexity}
                        """

## --------------------------------------------------------------------------------------##
## STAR mapping 
## --------------------------------------------------------------------------------------##
## Genome mapping with STAR
	## Note: If you want to use dm6_5kbrepeat you need to make a new StarIndex for that genome & add it to the config file (this is NOT the same as the bowtie2 genome!)

rule align_all_PairedEnd:
	input:
		readOne = 'Fastq/{sample}_R1_trim.fastq.gz',
		readTwo = 'Fastq/{sample}_R2_trim.fastq.gz',
		idx = 'Star_Index/genomeParameters.txt'
	output:
		STAR = 'Bam/{sample}_Aligned.sortedByCoord.out.bam',				
	params:
		starVer = config['module']['starVer'],
		refGenomePath = config['index']['dm6']['STAR'],
		sampleStr = "Bam/{sample}_",
	threads: 20
	shell:
		"""
		module purge && module load {params.starVer}
		STAR --genomeDir {params.refGenomePath} --runThreadN 20  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.sampleStr} --readFilesIn {input.readOne} {input.readTwo}
		"""

rule IndexBam:
	input:
		'Bam/{sample}_Aligned.sortedByCoord.out.bam'
	output:
		idx = 'Bam/{sample}_Aligned.sortedByCoord.out.bam.bai'
	params:
		samtoolsVer = config['module']['samtoolsVer']
	threads: 4
	shell:
		"""
		module purge && module load {params.samtoolsVer}
		samtools index {input}
		"""

rule bamToBed: 
	input:
		bam = "Bam/{sample}_Aligned.sortedByCoord.out.bam",
		idx = "Bam/{sample}_Aligned.sortedByCoord.out.bam.bai"
	output: 
		"Bed/{sample}_Aligned.sortedByCoord.out.bed"
	params:
		bedtoolsVer = config['module']['bedtoolsVer']
	shell:
		"""
		module purge && module load {params.bedtoolsVer}
		bedtools bamtobed -bedpe -i {input.bam} > {output} 
		"""

## Convert BAM files to bigWig
        ## parameters will need to be changed based on using a stranded/reverse stranded RNA library prep kit
rule BamtoBigWig:
		input:
				bam = 'Bam/{sample}_Aligned.sortedByCoord.out.bam',
				idx = 'Bam/{sample}_Aligned.sortedByCoord.out.bam.bai'
		output:
				bigWig = 'BigWig/{sample}_Aligned.sortedByCoord.bw',
				bigWigFwd = 'BigWig/{sample}_Aligned_fwd.sortedByCoord.bw',
				bigWigRev = 'BigWig/{sample}_Aligned_rev.sortedByCoord.bw'
		params:
				module = config['module']['deeptoolsVer'],
				genomeSize = config['general']['genomeSize'],
		threads: 4
		shell:
				"""
				module purge && module load {params.module}
				bamCoverage -split -i {input.bam} -p {threads}  --outFileFormat bigwig  -o {output.bigWig} --effectiveGenomeSize {params.genomeSize}
				bamCoverage -split -i {input.bam} -p {threads}  --outFileFormat bigwig  -o {output.bigWigRev} --effectiveGenomeSize {params.genomeSize} --filterRNAstrand reverse
				bamCoverage -split -i {input.bam} -p {threads}  --outFileFormat bigwig  -o {output.bigWigFwd} --effectiveGenomeSize {params.genomeSize} --filterRNAstrand forward
				"""

## --------------------------------------------------------------------------------------##
## Read Counting
## --------------------------------------------------------------------------------------##
## Read summarization using featureCounts
rule featureCounts:
	input:
		expand('Bam/{sample}_Aligned.sortedByCoord.out.bam', sample = sampleList)
	output:
		'featureCounts/counts_Aligned.txt'
	params:
		module = config['module']['subreadVer'],
		annotation = config['general']['annotation']
	threads: 4
	shell:
		"""
		module purge && module load {params.module}
		featureCounts -p -a {params.annotation} -o {output} -T 4 -t exon -g gene_id {input}
		"""

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")

## --------------------------------------------------------------------------------------##
## Data Visualization
## --------------------------------------------------------------------------------------##

## Use Bed and BigWig files for data visualization in IGV
## Data analysis done with DESeq2 in R. Use an OnDemand R Studio Interactive Session. 
## Will use the output featureCounts .txt file as an input for R analysis
## Use sampleSheet as metadata file for R analysis. 
