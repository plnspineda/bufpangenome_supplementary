#!/bin/bash
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=500GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

date

module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load BCFtools/1.17-GCC-11.2.0

files=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23/vg_deconstruct/concat_split_chr_1/*.vg.bcfmerge.vcf.gz)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo "Processing $file ..."

base=$(basename "$file" .vcf.gz)
echo "$base"

echo "Step: vcfbub"
mkdir -p clean_vcfbub || exit
##vcfbub 0.1.0
/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vcfbub/vcfbub -l 0 -r 1000000 --input "$file" | bgzip -@ 8 > clean_vcfbub/"$base".r1M.vcfbub.vcf.gz

echo "Step: vcfwave"
mkdir -p clean_vcfwave || exit
tabix -f -p vcf clean_vcfbub/"$base".r1M.vcfbub.vcf.gz
/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vcflib-1.0.12/vcflib/build/vcfwave -t 8 -L 1000000 -I 1000 clean_vcfbub/"$base".r1M.vcfbub.vcf.gz > clean_vcfwave/"$base".r1M.vcfbub.vcfwave.vcf

date
## next step is to normalise using bcftools