#!/usr/bin/env python3
"""
Filter MoChA calls to high-quality CNVs and flag event extent, using
MoChA's own `type` and arm-status (`p_arm`/`q_arm`) columns directly
instead of re-deriving classification from bdev/rel_cov.

Assumes df is a polars DataFrame from {pfx}.mca.calls.tsv with columns:
    sample_id, chrom, type, lod_baf_phase, p_arm, q_arm, rel_cov, length
(drop/rename the length filter block if you don't have that column)

`type` values follow MoChA convention: "Loss", "Gain", "CN-LOH",
"Undetermined", "CNP_Loss", "CNP_Gain" (germline CNP regions).
"""

import polars as pl
from mocha_calls_schema import read_mca_calls

ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def flag_mca_calls(
    df: pl.DataFrame,
    keep_types: tuple[str, ...] = ("Loss", "Gain"),
    drop_small_germline_dups: bool = True,
    small_dup_len_bp: float = 500_000,
    small_dup_rel_cov: float = 2.5,
) -> pl.DataFrame:
    chrom_str = pl.col("chrom").cast(pl.Utf8)
    is_acro = chrom_str.is_in(ACROCENTRIC)

    # Whole-chromosome event: both arms telomere-reaching, or for
    # acrocentric chromosomes, p_arm=="C" (centromere-reaching) proxy.
    whole_chrom = (
        ((pl.col("p_arm") == "T") & (pl.col("q_arm") == "T"))
        | (is_acro & (pl.col("p_arm") == "C") & (pl.col("q_arm") == "T"))
    )
    p_arm_event = (pl.col("p_arm") == "T") & (pl.col("q_arm") != "T") & ~whole_chrom
    q_arm_event = (pl.col("q_arm") == "T") & (pl.col("p_arm") != "T") & ~whole_chrom

    out = df.filter(
        pl.col("type").is_in(list(keep_types))     # keeps Loss/Gain, drops CN-LOH,   
    )

    if drop_small_germline_dups and "length" in df.columns:
        # MoChA's recommended germline-duplication filter (README, "Filter
        # callset"): a call passes if rel_cov is unremarkable (<2.1), OR
        # bdev is very low regardless of size, OR bdev/rel_cov are both
        # modest for a mid-size call (>500kb), OR bdev is a bit higher but
        # the call is large enough (>5Mb) that a germline dup this size
        # with this much imbalance would be unusual.

        length = pl.col("length")
        bdev = pl.col("bdev")
        rel_cov = pl.col("rel_cov")
        germline_dup_ok = (
            (rel_cov < 2.1) | (bdev < 0.05)
            | ((length > 5e5) & (bdev < 0.1) & (rel_cov < 2.5))
            | ((length > 5e6) & (bdev < 0.15))
        )
        out = out.filter(germline_dup_ok)


    out = out.with_columns([
        whole_chrom.alias("Whole_Chromosome"),
        p_arm_event.alias("P_Arm_Event"),
        q_arm_event.alias("Q_Arm_Event"),
    ])
    return out

def format_out(df: pl.DataFrame) -> pl.DataFrame:
    new_df = df.rename({"sample_id" : "SampleID",
                        "chrom": "Chr", 
                        "beg_GRCh38": "Start",
                        "end_GRCh38": "End",
                        "length": "Length",
                        "type" : "Type",
                        "cf": "Cell_Fraction",
                        "rel_cov" : "Relative_Copy_Number",
                        "bdev": "BAF_Deviation",
                        "lod_baf_phase": "LOD_BAF_Phase_Score" })
    
    new_df = new_df.select(["SampleID",
                            "Chr",
                            "Start",
                            "End",
                            "Length",
                            "Type", 
                            "Cell_Fraction", 
                            "Relative_Copy_Number", 
                            "BAF_Deviation", 
                            "LOD_BAF_Phase_Score",
                            "P_Arm_Event",
                            "Q_Arm_Event",
                            "Whole_Chromosome"])
    
    new_df = new_df.with_columns(pl.when(pl.col("Type") == "Loss")
                                   .then(pl.lit('DEL'))
                                   .when(pl.col("Type") == "Gain")
                                   .then(pl.lit("DUP")).alias("Type"))

    return(new_df)

if __name__ == "__main__":
    import sys

    calls = sys.argv[1]
    df_calls = read_mca_calls(calls)
    flagged = flag_mca_calls(df_calls)
    out = format_out(flagged)
    out.write_csv("mosaicDB.tsv", separator="\t")