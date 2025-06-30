#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH --mem=64GB

date

set -euo pipefail

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0
module load HTSlib/1.18-GCC-12.3.0

vcf="PCC_SB_variants/run2.PCC_SB.genomeWide.filtered.INDELs.vcf.gz"
dir="PCC_SB_variant"
name="sw17_"

mkdir -p "$dir"

base=$(basename "$vcf" .vcf.gz)

echo "$vcf"

echo "Subsampling"
bcftools stats "$vcf" > "$base".stats
bcftools view -S samples.txt "$vcf" -o "$dir/$name${base}.vcf.gz" -Oz
tabix -p vcf "$dir/$name${base}.vcf.gz"
bcftools stats "$dir/$name${base}.vcf.gz"

echo "Processing all filter"
bcftools view -f PASS "$dir/$name${base}.vcf.gz" -Oz -o "$dir/$name${base}.PASS.vcf.gz"
bcftools stats "$dir/$name${base}.PASS.vcf.gz" > "$dir/$name${base}.PASS.stats"

echo "Processing no homozygous ref but for all samples"
#num_samples=$(bcftools query -l "$vcf" | wc -l)
bcftools view -e "count(GT=\"0/0\") = 17" "$dir/$name${base}.PASS.vcf.gz" -Oz -o "$dir/$name${base}.PASS.nohomo.vcf.gz"
bcftools stats "$dir/$name${base}.PASS.nohomo.vcf.gz" > "$dir/$name${base}.PASS.nohomo.stats"
tabix -p vcf "$dir/$name${base}.PASS.nohomo.vcf.gz"

echo "Processing AUTOSOMES only"
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 "$dir/$name${base}.PASS.nohomo.vcf.gz" -Oz -o "$dir/$name${base}.PASS.nohomo.AUTO.vcf.gz"
tabix -p vcf "$dir/$name${base}.PASS.nohomo.AUTO.vcf.gz"

bcftools stats "$dir/$name${base}.PASS.nohomo.AUTO.vcf.gz" > "$dir/$name${base}.PASS.nohomo.AUTO.stats"
bcftools stats -s - "$dir/$name${base}.PASS.nohomo.AUTO.vcf.gz" > "$dir/$name${base}.PASS.nohomo.AUTO.samp.stats"

echo "Processing AUTOSOMES without chromosome 12"
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23 "$dir/$name${base}.PASS.nohomo.vcf.gz" -Oz -o "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.vcf.gz"
tabix -p vcf "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.vcf.gz"

bcftools stats "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.vcf.gz" > "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.stats"
bcftools stats -s - "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.vcf.gz" > "$dir/$name${base}.PASS.nohomo.AUTO.nochr12.samp.stats"

date
