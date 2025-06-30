# Identifying variants of the pangenome

Github page: https://github.com/vgteam/vg

`vg toolkit` is a multifunction tool for working with pangenome graphs.

1. First is to deconstruct the graph using `scripts/vg_decon.sh`. This will give a `.vcf.gz` file.
2. Next is to split the chromosome so I can rename the contigs easily. (They are in PanSN-Spec format, so have to convert it to chromosome and change the coordinates). This is done using `scripts/split_chrom.sh`
3. Rename the contigs and merge the files using `scripts/vg_rename.sh`
4. Filter and clean the vcf file using vcflib tools `scripts/vcfbub_wave.sh`
5. Normalise and get the variants using `scripts/bcfnorm.sh`
6. Adjust the coordinates for each chromosomes using `scripts/adjust_coor.sh`