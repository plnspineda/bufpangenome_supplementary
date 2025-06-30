#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=64GB
# #SBATCH --array=214-233
#SBATCH --ntasks-per-core=1

echo "Running task ${ARRAY_TASK_ID}"

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif"
samp="Sample${ARRAY_TASK_ID}"

pack=(/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/gam_files/"$samp"/chr_1_"$samp".pack)
gbz_dir="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/shortreads_variants/vg_giraffe_s100k_perchr/autoindex"

declare -A ref_paths
ref_paths["chr_2"]="RVUO#0#1#piece0 RVUO#0#1#piece1 RVUO#0#1#piece2 RVUO#0#1#piece3 RVUO#0#1#piece4 RVUO#0#1#piece5 RVUO#0#1#piece6 RVUO#0#1#piece7 RVUO#0#1#piece8 RVUO#0#1#piece9 RVUO#0#1#piece10 RVUO#0#1#piece11 RVUO#0#1#piece12"
ref_paths["chr_3"]="RVUO#0#2#piece0 RVUO#0#2#piece1 RVUO#0#2#piece2 RVUO#0#2#piece3 RVUO#0#2#piece4 RVUO#0#2#piece5 RVUO#0#2#piece6 RVUO#0#2#piece7 RVUO#0#2#piece8 RVUO#0#2#piece9 RVUO#0#2#piece10 RVUO#0#2#piece11 RVUO#0#2#piece12 RVUO#0#2#piece13 RVUO#0#2#piece14 RVUO#0#2#piece15 RVUO#0#2#piece16 RVUO#0#2#piece17 RVUO#0#2#piece18 RVUO#0#2#piece19 RVUO#0#2#piece20 RVUO#0#2#piece21 RVUO#0#2#piece22"
ref_paths["chr_4"]="RVUO#0#3#piece0 RVUO#0#3#piece1 RVUO#0#3#piece2 RVUO#0#3#piece3 RVUO#0#3#piece4 RVUO#0#3#piece5 RVUO#0#3#piece6 RVUO#0#3#piece7 RVUO#0#3#piece8 RVUO#0#3#piece9 RVUO#0#3#piece10 RVUO#0#3#piece11 RVUO#0#3#piece12 RVUO#0#3#piece13 RVUO#0#3#piece14 RVUO#0#3#piece15 RVUO#0#3#piece16 RVUO#0#3#piece17 RVUO#0#3#piece18 RVUO#0#3#piece19 RVUO#0#3#piece20 RVUO#0#3#piece21 RVUO#0#3#piece22 RVUO#0#3#piece23 RVUO#0#3#piece24 RVUO#0#3#piece25 RVUO#0#3#piece26 RVUO#0#3#piece27"
ref_paths["chr_1"]="RVUO#0#4#piece0 RVUO#0#4#piece1 RVUO#0#4#piece2 RVUO#0#4#piece3 RVUO#0#4#piece4 RVUO#0#4#piece5 RVUO#0#4#piece6 RVUO#0#4#piece7 RVUO#0#4#piece8 RVUO#0#4#piece9 RVUO#0#4#piece10 RVUO#0#4#piece11 RVUO#0#4#piece12 RVUO#0#4#piece13 RVUO#0#4#piece14 RVUO#0#4#piece15 RVUO#0#4#piece16 RVUO#0#4#piece17 RVUO#0#4#piece18 RVUO#0#4#piece19 RVUO#0#4#piece20 RVUO#0#4#piece21 RVUO#0#4#piece22 RVUO#0#4#piece23 RVUO#0#4#piece24 RVUO#0#4#piece25 RVUO#0#9#piece0 RVUO#0#9#piece1 RVUO#0#9#piece2 RVUO#0#9#piece3 RVUO#0#9#piece4 RVUO#0#9#piece5 RVUO#0#9#piece6 RVUO#0#9#piece7 RVUO#0#9#piece8 RVUO#0#9#piece9 RVUO#0#9#piece10 RVUO#0#9#piece11 RVUO#0#9#piece12 RVUO#0#9#piece13 RVUO#0#9#piece14 RVUO#0#9#piece15"
ref_paths["chr_5"]="RVUO#0#5#piece0 RVUO#0#5#piece1 RVUO#0#5#piece2 RVUO#0#5#piece3 RVUO#0#5#piece4 RVUO#0#5#piece5 RVUO#0#5#piece6 RVUO#0#5#piece7 RVUO#0#5#piece8 RVUO#0#5#piece9 RVUO#0#5#piece10 RVUO#0#5#piece11 RVUO#0#5#piece12 RVUO#0#5#piece13 RVUO#0#5#piece14 RVUO#0#5#piece15 RVUO#0#5#piece16 RVUO#0#5#piece17 RVUO#0#5#piece18 RVUO#0#5#piece19 RVUO#0#5#piece20 RVUO#0#5#piece21 RVUO#0#5#piece22 RVUO#0#5#piece23 RVUO#0#5#piece24 RVUO#0#5#piece25 RVUO#0#5#piece26"
ref_paths["chr_6"]="RVUO#0#6#piece0 RVUO#0#6#piece1 RVUO#0#6#piece2 RVUO#0#6#piece3 RVUO#0#6#piece4 RVUO#0#6#piece5 RVUO#0#6#piece6 RVUO#0#6#piece7 RVUO#0#6#piece8 RVUO#0#6#piece9 RVUO#0#6#piece10 RVUO#0#6#piece11 RVUO#0#6#piece12 RVUO#0#6#piece13 RVUO#0#6#piece14 RVUO#0#6#piece15 RVUO#0#6#piece16"
ref_paths["chr_7"]="RVUO#0#7#piece0 RVUO#0#7#piece1 RVUO#0#7#piece2 RVUO#0#7#piece3 RVUO#0#7#piece4 RVUO#0#7#piece5 RVUO#0#7#piece6 RVUO#0#7#piece7 RVUO#0#7#piece8 RVUO#0#7#piece9 RVUO#0#7#piece10 RVUO#0#7#piece11"
ref_paths["chr_8"]="RVUO#0#8#piece0 RVUO#0#8#piece1 RVUO#0#8#piece2 RVUO#0#8#piece3 RVUO#0#8#piece4 RVUO#0#8#piece5 RVUO#0#8#piece6 RVUO#0#8#piece7 RVUO#0#8#piece8 RVUO#0#8#piece9 RVUO#0#8#piece10"
ref_paths["chr_9"]="RVUO#0#10#piece0 RVUO#0#10#piece1 RVUO#0#10#piece2 RVUO#0#10#piece3 RVUO#0#10#piece4 RVUO#0#10#piece5 RVUO#0#10#piece6 RVUO#0#10#piece7 RVUO#0#10#piece8 RVUO#0#10#piece9 RVUO#0#10#piece10 RVUO#0#10#piece11 RVUO#0#10#piece12 RVUO#0#10#piece13 RVUO#0#10#piece14 RVUO#0#10#piece15 RVUO#0#10#piece16"
ref_paths["chr_10"]="RVUO#0#11#piece0 RVUO#0#11#piece1 RVUO#0#11#piece2 RVUO#0#11#piece3 RVUO#0#11#piece4 RVUO#0#11#piece5 RVUO#0#11#piece6 RVUO#0#11#piece7 RVUO#0#11#piece8 RVUO#0#11#piece9 RVUO#0#11#piece10 RVUO#0#11#piece11"
ref_paths["chr_11"]="RVUO#0#12#piece0 RVUO#0#12#piece1 RVUO#0#12#piece2 RVUO#0#12#piece3 RVUO#0#12#piece4 RVUO#0#12#piece5 RVUO#0#12#piece6"
ref_paths["chr_12"]="RVUO#0#13#piece0 RVUO#0#13#piece1 RVUO#0#13#piece2 RVUO#0#13#piece3 RVUO#0#13#piece4 RVUO#0#13#piece5 RVUO#0#13#piece6 RVUO#0#13#piece7 RVUO#0#13#piece8 RVUO#0#13#piece9 RVUO#0#13#piece10 RVUO#0#13#piece11 RVUO#0#13#piece12 RVUO#0#13#piece13 RVUO#0#13#piece14 RVUO#0#13#piece15 RVUO#0#13#piece16 RVUO#0#13#piece17 RVUO#0#13#piece18 RVUO#0#13#piece19 RVUO#0#13#piece20 RVUO#0#13#piece21 RVUO#0#13#piece22 RVUO#0#13#piece23 RVUO#0#13#piece24 RVUO#0#13#piece25 RVUO#0#13#piece26 RVUO#0#13#piece27 RVUO#0#13#piece28 RVUO#0#13#piece29"
ref_paths["chr_13"]="RVUO#0#14#piece0 RVUO#0#14#piece1 RVUO#0#14#piece2 RVUO#0#14#piece3 RVUO#0#14#piece4"
ref_paths["chr_14"]="RVUO#0#15#piece0 RVUO#0#15#piece1 RVUO#0#15#piece2 RVUO#0#15#piece3 RVUO#0#15#piece4 RVUO#0#15#piece5 RVUO#0#15#piece6 RVUO#0#15#piece7"
ref_paths["chr_15"]="RVUO#0#16#piece0 RVUO#0#16#piece1 RVUO#0#16#piece2 RVUO#0#16#piece3 RVUO#0#16#piece4 RVUO#0#16#piece5 RVUO#0#16#piece6 RVUO#0#16#piece7 RVUO#0#16#piece8 RVUO#0#16#piece9 RVUO#0#16#piece10 RVUO#0#16#piece11 RVUO#0#16#piece12 RVUO#0#16#piece13 RVUO#0#16#piece14 RVUO#0#16#piece15 RVUO#0#16#piece16 RVUO#0#16#piece17 RVUO#0#16#piece18 RVUO#0#16#piece19 RVUO#0#16#piece20"
ref_paths["chr_16"]="RVUO#0#17#piece0 RVUO#0#17#piece1 RVUO#0#17#piece2 RVUO#0#17#piece3 RVUO#0#17#piece4 RVUO#0#17#piece5 RVUO#0#17#piece6 RVUO#0#17#piece7 RVUO#0#17#piece8 RVUO#0#17#piece9 RVUO#0#17#piece10 RVUO#0#17#piece11"
ref_paths["chr_17"]="RVUO#0#18#piece0 RVUO#0#18#piece1 RVUO#0#18#piece2 RVUO#0#18#piece3 RVUO#0#18#piece4 RVUO#0#18#piece5 RVUO#0#18#piece6 RVUO#0#18#piece7 RVUO#0#18#piece8 RVUO#0#18#piece9 RVUO#0#18#piece10 RVUO#0#18#piece11 RVUO#0#18#piece12 RVUO#0#18#piece13 RVUO#0#18#piece14 RVUO#0#18#piece15 RVUO#0#18#piece16 RVUO#0#18#piece17 RVUO#0#18#piece18"
ref_paths["chr_18"]="RVUO#0#19#piece0 RVUO#0#19#piece1 RVUO#0#19#piece2 RVUO#0#19#piece3 RVUO#0#19#piece4 RVUO#0#19#piece5"
ref_paths["chr_19"]="RVUO#0#20#piece0 RVUO#0#20#piece1 RVUO#0#20#piece2 RVUO#0#20#piece3 RVUO#0#20#piece4 RVUO#0#20#piece5 RVUO#0#20#piece6 RVUO#0#20#piece7 RVUO#0#20#piece8 RVUO#0#20#piece9 RVUO#0#20#piece10 RVUO#0#20#piece11 RVUO#0#20#piece12"
ref_paths["chr_20"]="RVUO#0#21#piece0 RVUO#0#21#piece1 RVUO#0#21#piece2 RVUO#0#21#piece3 RVUO#0#21#piece4"
ref_paths["chr_21"]="RVUO#0#22#piece0 RVUO#0#22#piece1 RVUO#0#22#piece2 RVUO#0#22#piece3 RVUO#0#22#piece4 RVUO#0#22#piece5"
ref_paths["chr_22"]="RVUO#0#23#piece0 RVUO#0#23#piece1 RVUO#0#23#piece2 RVUO#0#23#piece3 RVUO#0#23#piece4 RVUO#0#23#piece5"
ref_paths["chr_23"]="RVUO#0#24#piece0 RVUO#0#24#piece1"

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