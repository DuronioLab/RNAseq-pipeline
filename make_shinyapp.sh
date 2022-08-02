#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8g
#SBATCH -t 6:00:00

git clone https://github.com/DuronioLab/RNA_shinyapp.git

mkdir RNA_shinyapp/input_data
mkdir RNA_shinyapp/output_data

printf "sample\trep\tbaseName\tbam\n" > ./RNA_shinyapp/input_data/sampleSheet.tsv

fileList=($(ls Bam/*_Aligned.sortedByCoord.out.bam))
for file in ${fileList[@]}
do
	tmp=${file##*/}
	basename=${tmp%%_Aligned*}
	rep=${basename##*_rep}
	sample=${basename%%_rep*}

	## Disabled for now b/c bam files can sometimes be too big
	
	#cp ./Bam/${tmp} ./RNA_shinyapp/input_data/${basename}_Aligned.sortedByCoord.out.bam
	#cp ./Bam/${tmp}.bai ./RNA_shinyapp/input_data/${basename}_Aligned.sortedByCoord.out.bai

	printf "${sample}\trep${rep}\t${basename}\tBam/${basename}_Aligned.sortedByCoord.out.bam\n" >> ./RNA_shinyapp/input_data/sampleSheet.tsv
	
done

cp ./featureCounts/counts_Aligned.txt ./RNA_shinyapp/input_data/featureCounts.txt

cp /proj/droniolb/genomeFiles/dm6/dmel-all-r6.20.ucsc.gtf ./RNA_shinyapp/input_data/dmel-all-r6.20.ucsc.gtf
cp /proj/droniolb/genomeFiles/dm6/fb_synonym_fb_2021_03.tsv ./RNA_shinyapp/input_data/fb_synonym_fb_2021_03.tsv

tar -zcvf RNA_shinyapp.tar.gz RNA_shinyapp