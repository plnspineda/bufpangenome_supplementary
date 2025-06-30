# Annotation of the structural variants, SNPs and InDels

## Repeat Masker

1. Take the sequences from the vcf file using `scripts/convert_vcf_to_fa.sh`
2. Run repeatmasker using:

    lib="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/repeatmasker/RepeatMasker/custom_lib/RepeatMaskerLib.h5.bostaurus-rm.withsats.fa"

    echo "repeatmasker for BIGIndel"
    fa="clean_SV_50bp_filtered.fa"
    base=$(basename $fa .fa)
    /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/repeatmasker/RepeatMasker/RepeatMasker -pa 24 -gff -lib "$lib" -dir final_"$base"_repeatmasker $fa -e ncbi

    date

## Variant Effect Predictor

To annotate the variants using vep: `scripts/vep.sh`

Further analysis were written in Rscripts which can be found here [link](https://github.com/plnspineda/bufpangenome_supplementary/tree/main/Rscripts)