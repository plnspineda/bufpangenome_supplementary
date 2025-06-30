#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=32GB
#SBATCH --array=10
#SBATCH --ntasks-per-core=1

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif"

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/pggb_2025.02.25/all_s100k/chr*p95_s100k/*smooth.final.gfa)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "Processing $file..."

base=$(basename "$file" .f3ab435.11fba48.e859484.smooth.final.gfa)
tmp=tmpdir_"$base"

mkdir -p $tmp
singularity exec "$sif" vg autoindex --workflow giraffe -g "$file" -p "$base" -t 16 -T $tmp

date