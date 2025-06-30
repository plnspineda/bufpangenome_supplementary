#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --array=0
#SBATCH --ntasks-per-core=1

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif"
samp="Sample230"

gams=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/gam_files/"$samp"/chr_1_*gam)
gam=${gams[$SLURM_ARRAY_TASK_ID]}
base=$(basename "$gam" .gam)
dir=$(dirname "$gam")
chr=$(basename "$gam" _"$samp".gam)
TMPDIR="$dir"
echo "Processing $gam"

gbz_dir="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/autoindex"
gbz="$gbz_dir"/"$chr".fa.giraffe.gbz

echo "STEP: VG PACK"
singularity exec "$sif" vg pack -x "$gbz" -g "$gam" --min-mapq 10 --threads 8 -o "$dir"/"$base".pack

echo "STEP: VG CALL"
singularity exec "$sif" vg call "$gbz" -k "$dir"/"$base".pack --sample "$base" --genotype-snarls --all-snarls --threads 8 > "$dir"/"$base".vcf

#if [ -f "${base}.pack" ] && [ -s "${base}.pack" ]; then
#    if rm "$gam" 2>/dev/null; then
#        echo "$gam removed"
#    else
#        echo "Error: Failed to remove $gam"
#    fi
#else
#    echo "${base}.pack does not exist or is empty, $gam not removed"
#fi

date