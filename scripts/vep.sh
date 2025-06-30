#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=72:00:00
#SBATCH --mem=64GB

module purge
module use /apps/modules/all
module load Singularity/3.10.5

vep="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/ensembl-vep.sif"
seq="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/variant_effect_predictor/database/sequences.fa.gz"
gene="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/variant_effect_predictor/database/genes.gff.gz"
vcf="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/vg_decon_2025.03.23/cleanvcfwave_adjusted_vcf/isec_output_dir/SV_50bp_filtered.vcf"

base="$(basename $vcf .vcf.gz)"

singularity exec $vep vep --fork 12 -i "$vcf" -o vep_ncbi_SV_final/"$base".res.txt --gff $gene --fasta $seq