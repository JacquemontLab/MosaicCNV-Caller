#!/usr/bin/env bash
# make_sampleFile.sh
# Usage: ./make_sampleFile.sh /path/to/vcf_dir [output_sampleFile.csv]

set -euo pipefail

VCF_DIR="${1:?Usage: $0 <vcf_dir> [output_file]}"
OUT="${2:-sampleFile.csv}"

if [[ ! -d "$VCF_DIR" ]]; then
    echo "Error: directory '$VCF_DIR' not found" >&2
    exit 1
fi

echo -e "Path,Chr" > "$OUT"

shopt -s nullglob
found=0
for f in "$VCF_DIR"/*.vcf "$VCF_DIR"/*.vcf.gz; do
    fname=$(basename "$f")

    # Extract chromosome token after chr/chromosome (case-insensitive): digits or X/Y
    if [[ "$fname" =~ [Cc][Hh][Rr](omosome)?[_.]?0*([0-9]+|[XYxy]) ]]; then
        chrom="chr${BASH_REMATCH[2]^^}"   # force lowercase 'chr' + uppercase X/Y (numbers unaffected)
        # if you'd rather keep numbers as numbers with no case issue, this is already fine
    else
        echo "Warning: could not parse a chromosome from $fname, skipping" >&2
        continue
    fi

    echo -e "$(realpath "$f"),${chrom}" >> "$OUT"
    found=$((found+1))
done

if [[ $found -eq 0 ]]; then
    echo "No VCF files matched a chromosome pattern in $VCF_DIR" >&2
    exit 1
fi

# sort by chromosome (2nd column now), natural/version order: chr1..chr22, chrX, chrY
{
    head -n1 "$OUT"
    tail -n +2 "$OUT" | sort -t$'\t' -k2,2V
} > "${OUT}.sorted" && mv "${OUT}.sorted" "$OUT"

echo "Wrote $found entries to $OUT"