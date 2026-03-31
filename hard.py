import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def _validate_inputs(counts: pd.DataFrame, metadata: pd.DataFrame, design: str):
    """
    Check that inputs meet DESeq2 requirements.
    """
    if not isinstance(counts, pd.DataFrame) or not isinstance(metadata, pd.DataFrame):
        raise TypeError("Counts and metadata must be Pandas DataFrames.")
        
    if counts.empty or metadata.empty:
        raise ValueError("DataFrames cannot be empty.")
        
    if counts.isnull().any().any():
        raise ValueError("Counts matrix contains NaN values.")

    if not pd.api.types.is_numeric_dtype(counts.values.flatten()):
        try:
            counts = counts.apply(pd.to_numeric)
        except Exception:
            raise ValueError("Counts must be strictly numeric.")

    if (counts.values < 0).any():
        raise ValueError("Counts matrix contains negative numbers.")
        
    if not np.all(np.mod(counts.values, 1) == 0):
        counts = counts.astype(int)

    count_samples = list(map(str, counts.columns))
    meta_samples = list(map(str, metadata.index))
    if set(count_samples) != set(meta_samples):
        raise ValueError("Counts columns and metadata index must exactly match.")

    if not design.strip().startswith("~"):
        raise ValueError("Design formula must start with '~', e.g., '~ condition'")

    return counts, metadata, design

def deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str) -> pd.DataFrame:
    """
    Run DESeq2 on a gene count matrix and sample metadata.
    """
    counts, metadata, design = _validate_inputs(counts, metadata, design)

    # Point rpy2 to the local R library path
    ro.r('.libPaths(c("~/R/libs", .libPaths()))')

    base = importr("base")
    stats = importr("stats")
    deseq2_r = importr("DESeq2")

    with localconverter(ro.default_converter + pandas2ri.converter):
        counts_r = ro.conversion.py2rpy(counts)
        meta_r = ro.conversion.py2rpy(metadata)
        
        # Let R natively parse the design formula
        design_r = stats.as_formula(design)
        
        dds = deseq2_r.DESeqDataSetFromMatrix(
            countData=counts_r,
            colData=meta_r,
            design=design_r
        )
        dds = deseq2_r.DESeq(dds)
        res = deseq2_r.results(dds)
        
        return ro.conversion.rpy2py(base.as_data_frame(res))

if __name__ == "__main__":

    print("Running deseq2()...\n")
    counts_df = pd.read_json("dataset/airway_counts.json", orient="index")
    meta_df = pd.read_json("dataset/airway_metadata.json", orient="index")
    
    results = deseq2(counts_df, meta_df, "~ dex")
    print(results.head(5))
 