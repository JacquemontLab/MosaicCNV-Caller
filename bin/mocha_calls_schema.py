#!/usr/bin/env python3
"""
Schema and reader for MoChA's mosaic calls table ({pfx}.calls.tsv), per the
column descriptions in https://github.com/freeseek/mocha (README, "Call
chromosomal alterations" section).

MoChA writes missing values as the literal string "nan" rather than an
empty field, which polars' own float parser won't recognize by default
(and even where it would, we want it explicit and applied consistently
across every column, not just the numeric ones). read_mca_calls() below
loads everything as text first, swaps "nan" -> null, then casts.
"""

import re
import polars as pl

# Every column except beg_XXXXXX/end_XXXXXX, which carry the genome build
# in their suffix (e.g. beg_GRCh38, end_GRCh37) and so aren't fixed keys.
CALLS_SCHEMA: dict[str, pl.DataType] = {
    "sample_id": pl.Utf8,
    "computed_gender": pl.Utf8,     # M, F, U (unknown), or K (Klinefelter)
    "chrom": pl.Utf8,
    "length": pl.Int64,
    "p_arm": pl.Utf8,                # N = doesn't extend, C = to centromere, T = to telomere
    "q_arm": pl.Utf8,                # same coding as p_arm, for the long arm
    "n_sites": pl.Int64,
    "n_hets": pl.Int64,
    "n50_hets": pl.Int64,
    "bdev": pl.Float64,
    "bdev_se": pl.Float64,
    "rel_cov": pl.Float64,
    "rel_cov_se": pl.Float64,
    "lod_lrr_baf": pl.Float64,
    "lod_baf_phase": pl.Float64,
    "n_flips": pl.Int64,             # -1 sentinel: LRR+BAF model used instead of BAF+phase model
    "baf_conc": pl.Float64,
    "lod_baf_conc": pl.Float64,
    "type": pl.Utf8,                 # Loss, Gain, CN-LOH, Undetermined, or CNP_Loss / CNP_Gain
    "cf": pl.Float64,
}

CALLS_COLUMN_DESCRIPTIONS: dict[str, str] = {
    "sample_id": "Sample identifier.",
    "computed_gender": "Inferred sample sex/gender for this call: M, F, U (unknown), or K (Klinefelter).",
    "chrom": "Chromosome the call is on.",
    "beg_XXXXXX": "Call start position, in coordinates of genome build XXXXXX (e.g. beg_GRCh38).",
    "end_XXXXXX": "Call end position, in coordinates of genome build XXXXXX (e.g. end_GRCh38).",
    "length": "Call length in base pairs.",
    "p_arm": "Whether the call reaches into the short arm: N = no, C = reaches the centromere, T = reaches the telomere.",
    "q_arm": "Same coding as p_arm, for the long arm.",
    "n_sites": "Number of sites (variants) used to make the call.",
    "n_hets": "Number of heterozygous sites used to make the call.",
    "n50_hets": "N50 of the distances between consecutive heterozygous sites used in the call.",
    "bdev": "Estimated BAF deviation from 0.5 (the allelic-imbalance statistic).",
    "bdev_se": "Standard error of the bdev estimate.",
    "rel_cov": "Estimated relative coverage/copy number, derived from LRR or sequencing depth (diploid = 2).",
    "rel_cov_se": "Standard error of the rel_cov estimate.",
    "lod_lrr_baf": "LOD score for the call under the LRR+BAF model.",
    "lod_baf_phase": "LOD score for the call under the BAF+genotype-phase model.",
    "n_flips": "Number of phase flips detected for BAF+phase model calls; -1 if the LRR+BAF model was used instead.",
    "baf_conc": "BAF phase concordance across the call's phased heterozygous sites.",
    "lod_baf_conc": "Genome-wide-corrected LOD score for the BAF phase concordance model.",
    "type": "Call type inferred from LRR/relative coverage: Loss, Gain, CN-LOH, Undetermined, or CNP_Loss/CNP_Gain (germline copy-number polymorphism).",
    "cf": "Estimated cell fraction (from bdev+type, or rel_cov+type if bdev/bdev_se are missing).",
}

BEG_END_RE = re.compile(r"^(beg|end)_")


def read_mca_calls(path: str) -> pl.DataFrame:
    """
    Read a MoChA {pfx}.calls.tsv file, casting each column per
    CALLS_SCHEMA and converting the literal string "nan" to a real null
    in every column (not just numeric ones) before casting.
    """
    # Read as all-Utf8 first so "nan" tokens don't get silently mishandled
    # by type inference, and so the null-swap below applies uniformly.
    raw = pl.read_csv(path, separator="\t", infer_schema_length=0)

    exprs = []
    for col in raw.columns:
        dtype = pl.Int64 if BEG_END_RE.match(col) else CALLS_SCHEMA.get(col, pl.Utf8)
        expr = pl.when(pl.col(col) == "nan").then(None).otherwise(pl.col(col))
        if dtype != pl.Utf8:
            expr = expr.cast(dtype, strict=False)
        exprs.append(expr.alias(col))

    return raw.with_columns(exprs)


if __name__ == "__main__":
    import sys

    df = read_mca_calls(sys.argv[1])
    print(df.schema)
    print(df.null_count())