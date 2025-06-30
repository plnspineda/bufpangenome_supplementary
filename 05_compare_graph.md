# Compare linear vs graph

1. Extract specific paths by running vg call again using `scripts/vgcall_swpc.sh` for the Philippine swamp buffalo, and `scripts/vgcall_rvuo.sh` for river buffalo.
2. Compress the files:

        for vcf in Sample2*/*.vcf; do
        echo $vcf
        bgzip $vcf
        tabix -p vcf $vcf.gz
        done

3. Concatenate the files using `scripts/bcfconcat.sh`
4. Merge the files:

        for i in merged*.vcf.gz; do
            echo "$i"
            gunzip $i
            bgzip $(basename $i .gz)
            tabix -p vcf $i
        done

        ls merged*.vcf.gz > merged.list

        bcftools merge --file-list merged.list -Oz --threads 8 -o all_samples.merged.vcf.gz
        rm merged.list

5. Adjust the coordinates using `scripts/adjust_path.sh`
6. Get each variants from the `get_variant_scripts` folder
7. Compare linear vs graph using:

        v1="$1"
        v2="$2"
        dir="$3"

        bcftools isec "$v1" "$v2" -p isec_"$dir"
