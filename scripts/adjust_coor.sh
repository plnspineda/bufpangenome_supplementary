#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

date
set -euo pipefail

module purge
module use /apps/modules/all
module load HTSlib/1.18-GCC-12.3.0
module load Python/3.11.3-GCCcore-12.3.0
module load BCFtools/1.17-GCC-11.2.0

dir="cleanvcfwave_adjusted_vcf"
coord="swpc_coord.txt"
TMPDIR="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23"

mkdir -p "$dir"

for vcf in clean_vcfwave/*vcfbub.vcfwave.vcf; do
  echo "$vcf"
  base=$(basename "$vcf" .vcf)
  python adjust.py "$vcf" "$dir"/"$base".vcf "$coord" false
  bgzip "$dir"/"$base".vcf
  bcftools sort "$dir"/"$base".vcf.gz -Oz -o "$dir"/"$base".sorted.vcf.gz
  tabix -p vcf "$dir"/"$base".sorted.vcf.gz
done
date