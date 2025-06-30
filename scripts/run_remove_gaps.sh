#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=01:00:00
#SBATCH --mem=16GB

module use /apps/modules/all
module load Python/3.11.3-GCCcore-12.3.0
module load seqtk/1.3-GCC-11.2.0

for i in *fna; do
    echo "python 01_remove_gaps.py $i"
    python 01_remove_gaps.py "$i"
done