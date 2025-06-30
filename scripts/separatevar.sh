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

echo "Processing"
vcf="merged.18samples.vcf.gz"
base=$(basename "$vcf" .vcf.gz)

mkdir -p tmpdir

echo "Separate variants"
echo "big indels"
bcftools view -v indels "$vcf" | awk '($0 ~ /^#/) || (length($5) - length($4) > 50 || length($4) - length($5) > 50)' | bcftools view -H > "$base".BIGindel.vcf

echo "small indels"
bcftools view -v indels "$vcf" | awk '($0 ~ /^#/) || (length($5) - length($4) < 50 && length($4) - length($5) < 50)' | bcftools view -H > "$base".smallIndel.vcf

echo "SNPs"
bcftools view -v snps "$vcf" | bcftools view -H > "$base".SNPs.vcf

echo "SVs"
gunzip -c "$vcf" | \
awk 'BEGIN{OFS="\t"}
     /^#/ {
       print > "header.tmp";
       next
     }
     {
       ref_len = length($4)
       alt_len = length($5)

       # Skip SNVs and INDELs
       if ((ref_len == 1 && alt_len == 1) || (ref_len == 1 || alt_len == 1))
         next

       len_diff = (ref_len > alt_len ? ref_len - alt_len : alt_len - ref_len)

       if (len_diff > 50)
         print >> "sv_gt50bp.body"
       else if (len_diff < 50)
         print >> "sv_lt50bp.body"
     }

     END {
       system("cat header.tmp sv_gt50bp.body > merged.18samples.bigSV.vcf")
       system("cat header.tmp sv_lt50bp.body > merged.18samples.smallSV.vcf")
       system("rm header.tmp sv_gt50bp.body sv_lt50bp.body")
     }'

date
