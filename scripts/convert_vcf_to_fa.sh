#!/bin/bash

vcf_file="$1"
base=$(basename "$vcf_file" .vcf)
dir=$(dirname "$1")
output_fasta="$dir"/"$base".fa

touch "$output_fasta"

# Process the VCF
awk -F'\t' 'BEGIN { OFS="\t" }
  !/^#/ {
    chrom = $1
    pos = $2
    id = $3
    gsub(/>/, "_", id)
    ref = $4
    alt = $5

    # Output ref
    print ">" chrom "-" pos "-" id "-ref-len" length(ref) >> "'"$output_fasta"'"
    print ref >> "'"$output_fasta"'"

    # Output each alt
    split(alt, alt_alleles, ",")
    for (a in alt_alleles) {
      alt_seq = alt_alleles[a]
      print ">" chrom "-" pos "-" id "-alt" a "-len" length(alt_seq) >> "'"$output_fasta"'"
      print alt_seq >> "'"$output_fasta"'"
    }
  }
' "$vcf_file"