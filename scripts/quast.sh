#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

set -euo pipefail

module purge
module use /apps/modules/all
module load quast/4.5-foss-2016uofa-Python-2.7.11

for i in genomes/*_ungapped.fa; do
  echo "$i"
  base=$(basename "$i" .fa)
  quast.py "$i" -o quast_"${base}"
done