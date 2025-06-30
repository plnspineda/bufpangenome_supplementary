#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

set -euo pipefail

date
module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load SAMtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/input-files/input_pggb/*.fa)
asm=${files[$SLURM_ARRAY_TASK_ID]}

base=$(basename "$asm" .fa)

mkdir -p $base || cd $base
#samtools faidx "$asm"

singularity run /hpcfs/users/a1812753/buffalo_pangenome/pggb_latest.sif pggb -i "$asm" -n 11 -p 95 -s 100k -o "$base"_p95_s100k

date