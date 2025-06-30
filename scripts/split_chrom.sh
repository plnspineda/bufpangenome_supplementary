#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=12GB
#SBATCH --array=10
#SBATCH --ntasks-per-core=1

date

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23/vg_deconstruct/*.vcf.gz)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "$file"

dir=$(dirname "$file")
chr=$(basename "$file" | cut -d'.' -f1)

mkdir -p "$dir"/split_"$chr" || exit

dir=$(dirname "$file")
chr=$(basename "$file" | cut -d'.' -f1)
echo "$dir"
echo "$chr"

mkdir -p "$dir/split_$chr" || exit

for sample in $(bcftools query -l "$file"); do
    base_name=$(basename "$file" .vcf.gz)
    bcftools view --threads 8 -c1 -Oz -s "$sample" -o "$dir/split_$chr/${base_name}.$sample.vcf.gz" "$file"
done

echo "Done!"
date