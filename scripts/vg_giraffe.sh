#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --array=10
#SBATCH --ntasks-per-core=1

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif"

samp="Sample230"

gbz_array=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/autoindex/chr_*.fa.giraffe.gbz)
gbz=${gbz_array[$SLURM_ARRAY_TASK_ID]}
echo "Processing chromosome $gbz..."
gbzdir=$(dirname "$gbz")

export TMPDIR=$PWD
echo "temporary directory: $TMPDIR"

sample=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads/ph_samples/CleanedReads/"$samp".R1_val_1.fq.gz
echo "Processing shortreads $sample..."

sample_base=$(basename "$sample" .R1_val_1.fq.gz)
sample_dir=$(dirname "$sample")
r1="$sample_dir"/"$sample_base".R1_val_1.fq.gz
r2="$sample_dir"/"$sample_base".R2_val_2.fq.gz

gbzbase=$(basename "$gbz" .giraffe.gbz)
base=$(basename "$gbz" .fa.giraffe.gbz)

echo "read1: $r1"
echo "read2: $r2"
echo "files: $gbzbase"

mkdir -p $sample_base

#map shortreads
singularity exec "$sif" vg giraffe \
  -p \
  -t 8 \
  -Z "$gbzdir"/"$gbzbase"*.giraffe.gbz \
  -d "$gbzdir"/"$gbzbase"*.dist \
  -m "$gbzdir"/"$gbzbase"*.min \
  -f "$r1" -f "$r2" > $sample_base/"$base"_"$sample_base".gam

singularity exec "$sif" vg stats -a $sample_base/"$base"_"$sample_base".gam
date