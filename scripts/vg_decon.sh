#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

date

module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load BCFtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/pggb_2025.02.25/all_s100k/chr_1_p95_s100k/*.gfa)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "Processing $file ..."

base=$(basename "$file" .gfa)
echo "gfa file path: $file"
echo "output file: $base.vcf.gz"

echo "Step1: vg deconstruct"
mkdir -p vg_deconstruct || exit
##vg_v1.64.0
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg deconstruct --threads 8 -P SWPC --all-snarls "$file" | bgzip -c -@ 8 > vg_deconstruct/"$base".vcf.gz

date