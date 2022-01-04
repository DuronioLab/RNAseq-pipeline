#pipepath=rna-seq-pipeline

snakeFile=$1

if [[ ! -d ./Err_Out/ ]]; then
	mkdir ./Err_Out/
fi

if [[ -z ${snakeFile} ]]; then
	snakeFile='./snakefile_RNA_snakemake.py'
else
	snakeFile=${snakeFile}
fi

echo ${snakeFile}

snakemake --printshellcmds -nrs --snakefile ${snakeFile}

#snakemake --rerun-incomplete --printshellcmds -s ${snakeFile} --latency-wait 60 --cluster "sbatch -J {rule} -o ./Err_Out/slurm-%j.out -e ./Err_Out/slurm-%j.err -N1 --time 2:00:00 --mem=16G" --jobs 500
snakemake --cluster-config slurmConfig --printshellcmds -s ${snakeFile} --rerun-incomplete --latency-wait 100000 --cluster  "sbatch -J {rule} -o ./Err_Out/slurm-%j.out -e ./Err_Out/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" --jobs 500
