#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --array=223,231,232
#SBATCH --ntasks-per-core=1

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

samp="Sample${SLURM_ARRAY_TASK_ID}"
echo "Processing $samp"

base=$(basename "$samp")

ls "$samp"/*.vcf.gz > tmp_"$base".list
bcftools concat --file-list tmp_"$base".list -Oz --threads 8 -o merged."$base".vcf.gz
rm tmp_"$base".list