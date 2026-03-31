import pytest
import pandas as pd
from hard import deseq2

@pytest.fixture
def dummy_data():
    """
    Create a robust dummy dataset for testing the DESeq2 wrapper.
    """
    counts = pd.DataFrame({
        "sample1": [100, 200, 50, 10],
        "sample2": [120, 190, 45, 12],
        "sample3": [400, 50, 500, 100],
        "sample4": [450, 60, 480, 110],
        "sample5": [110, 220, 48, 15],
        "sample6": [380, 55, 520, 95]
    }, index=["gene1", "gene2", "gene3", "gene4"])
    
    metadata = pd.DataFrame({
        "condition": ["control", "control", "treatment", "treatment", "control", "treatment"]
    }, index=["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])
    
    return counts, metadata

def test_deseq2_runs_successfully(dummy_data):
    counts, metadata = dummy_data
    res = deseq2(counts, metadata, "~ condition")
    
    assert isinstance(res, pd.DataFrame)
    assert not res.empty
    
    expected_cols = {"baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"}
    assert expected_cols.issubset(res.columns)
    assert list(res.index) == list(counts.index)

def test_input_validation_empty():
    counts = pd.DataFrame()
    metadata = pd.DataFrame()
    
    with pytest.raises(ValueError, match="DataFrames cannot be empty."):
        deseq2(counts, pd.DataFrame({"cond":["A"]}), "~ cond")
        
    with pytest.raises(ValueError, match="DataFrames cannot be empty."):
        deseq2(pd.DataFrame({"S1":[1]}), metadata, "~ cond")

def test_input_validation_matching_names(dummy_data):
    counts, metadata = dummy_data
    metadata.index = ["S1", "S2", "S3", "S4", "S5", "S6"]
    
    with pytest.raises(ValueError, match="Counts columns and metadata index must exactly match."):
        deseq2(counts, metadata, "~ condition")

def test_input_validation_non_negative(dummy_data):
    counts, metadata = dummy_data
    counts.loc["gene1", "sample1"] = -5
    
    with pytest.raises(ValueError, match="Counts matrix contains negative numbers."):
        deseq2(counts, metadata, "~ condition")
