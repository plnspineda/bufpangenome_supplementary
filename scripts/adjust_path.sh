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

dir="adjusted_coord_vcf"
coord="swpc_coord.txt"
TMPDIR="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/swamp_asref_vcf_2"
adjust="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/scripts/adjust_swpc.py"

vcf="all_samples.merged.vcf.gz"
echo "$vcf"
base=$(basename "$vcf" .vcf.gz)
python "$adjust" "$vcf" "$dir"/"$base".vcf "$coord" false
bgzip "$dir"/"$base".vcf
bcftools sort "$dir"/"$base".vcf.gz -Oz -o "$dir"/"$base".sorted.vcf.gz
tabix -p vcf "$dir"/"$base".sorted.vcf.gz

date