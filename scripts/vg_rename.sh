#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mem=12GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

date

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23/vg_deconstruct/split_chr_1)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "$file"

if [[ -z "$file" || ! -d "$file" ]]; then
  echo "Error: '$file' is not a valid directory" >&2
  exit 1
fi

basedir=$(basename "$file")
dir=$(dirname "$file")/renamed_"$basedir"

for vcf in "$file"/*.vcf.gz; do
  echo "$vcf"
  base=$(basename "$vcf" .vcf.gz).renamed.vcf.gz
  mkdir -p "$dir"
  zcat "$vcf" | awk 'BEGIN {OFS="\t"} /^#CHROM/ {for (i=1; i<=NF; i++) gsub(/#[0-9]+#piece[0-9]+/, "", $i)} {print}' | bgzip > "$dir"/"$base"
done

asm=(RVAZ#1 RVAZ#2 RVNR#1 RVNR#2 RVUO#0 RVND#1 RVND#2 RVCU#0 SWPC#0 SWCU#0 SWWA#0)

if [[ ! -d "$dir" ]]; then
  echo "Error: Directory '$dir' does not exist" >&2
  exit 1
fi

for line in "${asm[@]}"; do
  echo "$line"
  ls "$dir"/*"$line"*.renamed.vcf.gz > "$dir"/"$line".list

  if [[ -s "$dir"/"$line".list ]]; then
    concatdir=$(dirname "$file")/concat_"$basedir"
    mkdir -p "$concatdir"
    bcftools concat -f "$dir"/"$line".list --naive -o "$concatdir"/"$line".concat.vcf.gz -O z
    bcftools sort "$concatdir"/"$line".concat.vcf.gz -Oz > "$concatdir"/"$line".sort.concat.vcf.gz
    tabix -p vcf "$concatdir"/"$line".sort.concat.vcf.gz
  else
    echo "Warning: No files found for $line. Skipping bcftools concat."
  fi
done


ls "$concatdir"/*.sort.concat.vcf.gz > "$concatdir"/merge.list
chr="${basedir#split_}"

bcftools merge --threads 4 -l "$concatdir"/merge.list -Oz -o "$concatdir"/"$chr".vg.bcfmerge.vcf.gz

echo "Done!"
date