# Preparation of the genome assemblies input files

1. First is to determine the homologous chromosomes using `scripts/mashmap.sh`
2. Determine contig N50 using `scripts/quast.sh`
3. Next is to separate each chromosomes `scripts/separate_chr_contig.sh`
4. Then remove the gaps `scripts/run_remove_gaps.sh`
5. Then concatenate the contigs per chromosome in all assemblies `scripts/concat_chr.sh`
6. I have to manually do it for chromosome 1 with

        cat RVCU04_ungapped.fa RVCU09_ungapped.fa RVND24_ungapped.fa RVND29_ungapped.fa RVND14_ungapped.fa RVND19_ungapped.fa RVUO04_ungapped.fa RVUO09_ungapped.fa RVAZ14_ungapped.fa RVAZ19_ungapped.fa RVAZ24_ungapped.fa RVAZ29_ungapped.fa RVNR14_ungapped.fa RVNR19_ungapped.fa RVNR24_ungapped.fa RVNR29_ungapped.fa SWCU01_ungapped.fa SWPC01_ungapped.fa SWWA01_ungapped.fa > chr_1.fa

7. The input files for PGGB would be per chromosome containing sequences from each assemblies with headers like `RVCU#0#4#piece0` from PanSN-spec (https://github.com/pangenome/PanSN-spec)

## Scaffolding the Pakistan river assemblies

To scaffold the Pakistan river assemblies, I use `scripts/ragtag.sh`. Then run `scripts/dotplot.sh` to visualise the scaffolding using dotplot.
