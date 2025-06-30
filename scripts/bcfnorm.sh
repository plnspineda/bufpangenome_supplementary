#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

date

module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load BCFtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23/clean_vcfwave_chr4/chr_4.*r1M.vcfbub.vcfwave.vcf)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "Processing $file ..."

base=$(basename "$file" .vcf.gz)
echo "vcf file path: $file"

## get SNPs
echo "Step: SNPs"
mkdir -p tmpdir || exit
mkdir -p vcf_SNPs
bcftools view -m2 -M2 -v snps "$file" | \
bcftools norm --rm-dup exact | \
bcftools norm -m -any | \
bcftools sort --temp-dir tmpdir -o vcf_SNPs/"$base".SNP.vcf.gz --output-type z

## get indels
echo "Step: indels"
mkdir -p vcf_indels
bcftools view -v indels -i 'abs(ILEN) > 0 && abs(ILEN) < 50' "$file" | \
bcftools norm --rm-dup exact | \
bcftools norm -m -any | \
bcftools sort --temp-dir tmpdir -o vcf_indels/"$base".indel.vcf.gz --output-type z

## get big indels
echo "Step: big indels"
mkdir -p vcf_big_indels
bcftools view -v indels -i 'abs(ILEN) > 50' "$file" | \
bcftools norm --rm-dup exact | \
bcftools norm -m -any | \
bcftools sort --temp-dir tmpdir -o vcf_big_indels/"$base".BIGindel.vcf.gz --output-type z

## get small SV
echo "Step: small SV"
mkdir -p vcf_smallSV
bcftools view -v other -i 'abs(ILEN) > 0 && abs(ILEN) < 50' "$file" | \
bcftools norm --rm-dup exact | \
bcftools norm -m -any | \
bcftools sort --temp-dir tmpdir -o vcf_smallSV/"$base".smallSV.vcf.gz --output-type z

## get SV
echo "Step: SV"
mkdir -p vcf_SV
bcftools view -v other -i 'abs(ILEN) > 50' "$file" | \
bcftools norm --rm-dup exact | \
bcftools norm -m -any | \
bcftools sort --temp-dir tmpdir -o vcf_SV/"$base".SV.vcf.gz --output-type z

echo "Done!"
date