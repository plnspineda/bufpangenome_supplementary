#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=72:00:00
#SBATCH --mem=64GB
#SBATCH --array=223
#SBATCH --ntasks-per-core=1

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif"
samp="Sample${SLURM_ARRAY_TASK_ID}"

# Define chromosome to ref-path mapping based on .fai data
declare -A ref_paths
ref_paths["chr_1"]="SWPC#0#1#piece0 SWPC#0#1#piece1 SWPC#0#1#piece2 SWPC#0#1#piece3"
ref_paths["chr_2"]="SWPC#0#2#piece0 SWPC#0#2#piece1"
ref_paths["chr_3"]="SWPC#0#3#piece0 SWPC#0#3#piece1"
ref_paths["chr_4"]="SWPC#0#4#piece0 SWPC#0#4#piece1 SWPC#0#4#piece2 SWPC#0#4#piece3 SWPC#0#4#piece4 SWPC#0#4#piece5"
ref_paths["chr_5"]="SWPC#0#5#piece0"
ref_paths["chr_6"]="SWPC#0#6#piece0 SWPC#0#6#piece1 SWPC#0#6#piece2"
ref_paths["chr_7"]="SWPC#0#7#piece0"
ref_paths["chr_8"]="SWPC#0#8#piece0"
ref_paths["chr_9"]="SWPC#0#9#piece0 SWPC#0#9#piece1"
ref_paths["chr_10"]="SWPC#0#10#piece0"
ref_paths["chr_11"]="SWPC#0#11#piece0"
ref_paths["chr_12"]="SWPC#0#12#piece0"
ref_paths["chr_13"]="SWPC#0#13#piece0"
ref_paths["chr_14"]="SWPC#0#14#piece0"
ref_paths["chr_15"]="SWPC#0#15#piece0"
ref_paths["chr_16"]="SWPC#0#16#piece0"
ref_paths["chr_17"]="SWPC#0#17#piece0 SWPC#0#17#piece1"
ref_paths["chr_18"]="SWPC#0#18#piece0"
ref_paths["chr_19"]="SWPC#0#19#piece0"
ref_paths["chr_20"]="SWPC#0#20#piece0"
ref_paths["chr_21"]="SWPC#0#21#piece0"
ref_paths["chr_22"]="SWPC#0#22#piece0"
ref_paths["chr_23"]="SWPC#0#23#piece0 SWPC#0#23#piece1"

pack=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/gam_files/"$samp"/chr_*_"$samp".pack)
gbz_dir="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/autoindex"

if [ ${#pack[@]} -eq 0 ]; then
    echo "Error: No .pack files found for $samp"
    exit 1
fi

mkdir -p "$samp"
dir="$samp"

for pack in "${pack[@]}"; do
    base=$(basename "$pack" .pack)
    chr=$(basename "$pack" _"$samp".pack)
    TMPDIR="$dir"
    echo "Processing $pack for $samp (chromosome: $chr)"

    gbz="$gbz_dir"/"$chr".fa.giraffe.gbz
    if [ ! -f "$gbz" ]; then
        echo "Error: GBZ file $gbz does not exist"
        exit 1
    fi

    ref_path_args=()
    if [ -n "${ref_paths[$chr]}" ]; then
        for path in ${ref_paths[$chr]}; do
            ref_path_args+=("-p" "$path")
        done
    else
        echo "Error: No ref-paths defined for $chr"
        exit 1
    fi

    echo "STEP: VG CALL with ref-paths: ${ref_path_args[*]}"
    echo "singularity exec $sif vg call $gbz -k $pack ${ref_path_args[*]} --sample $samp --genotype-snarls --all-snarls --threads 8 > $dir/$base.vcf"
    singularity exec "$sif" vg call "$gbz" -k "$pack" "${ref_path_args[@]}" --sample "$samp" --genotype-snarls --all-snarls --threads 8 > "$dir"/"$base".vcf
done

date