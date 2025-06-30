# Mapping short reads with pangenome graph

1. Index the graph using `scripts/vg_autoindex.sh`
2. align the cleaned shortreads using `scripts/vg_giraffe.sh`
3. compact the aligns and call the variants using `scripts/vg_call.sh`
4. Merge and concatenate the vcf files
    - compress the vcf files

        for vcf in Sample2*/*.vcf; do
            echo "$vcf"
            bgzip "$vcf"
            tabix -p vcf "${vcf}.gz"
        done
    
    - concatenate

        for vcf in Sample*/chr_*.vcf.gz; do
            [ -f "$vcf" ] || continue
            out="$(dirname "$vcf")/renamed_$(basename "$vcf" .vcf.gz).vcf.gz"
            bcftools reheader --samples rename_list.txt -o "$out" "$vcf"
            tabix -p vcf -f "$out"
        done

        for samp in Sample*/; do
            echo $samp
            ls "$samp"renamed_*.vcf.gz > "$samp"tmp.txt
            bcftools concat --file-list "$samp"tmp.txt --threads 4 -O z -o "$samp"/merged.vcf.gz
            rm -f "$samp"tmp.txt
        done

    - filter and normalise

        bcftools view -f PASS "$vcf" -Oz -o "$base".PASS.vcf.gz
        tabix -p vcf "$base".PASS.vcf.gz

        # normalise
        echo "Normalise vcf file"
        bcftools norm --rm-dup exact "$base".PASS.vcf.gz | \
        bcftools norm -m -any | \
        bcftools sort --temp-dir tmpdir -o "$base".NORM.vcf.gz --output-type z
        tabix -p "$base".NORM.vcf.gz

    - merge the vcf files

        bcftools merge --file-list passed_all_merge_list.txt -Oz -o merged.samples.vcf.gz
        bcftools stats merged.samples.vcf.gz > merged.samples.stats
        bcftools stats -s - merged.samples.vcf.gz > merged.samples.SAMP.stats

5. And then separate the variants using `scripts/separatevar.sh`