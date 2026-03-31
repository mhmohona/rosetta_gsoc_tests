# [Rosetta GSoC Code Tests](https://github.com/rstats-gsoc/gsoc2026/wiki/rosetta%3A-Python-wrappers-for-Bioconductor-packages-via-rpy2)

Solutions for the overarching GSoC 2026 project tests for the `rosetta` wrapper, which provides Python interfaces for Bioconductor differential expression engines via `rpy2`. The code prioritizes native R statistical proxy modeling rather than raw string execution, ensuring a clean and robust Python API that safely checks inputs before parsing.

## 1. Easy Task ([easy.R](https://github.com/mhmohona/rosetta_gsoc_tests/blob/main/easy.R))

Native R code that models the well-known `airway` dataset inside a `DESeqDataSet` object, factoring by cell structures and dex treatment (`~ cell + dex`). The output cleanly prints the top 5 differentially expressed genes calculated by lowest false discovery rate.

**Output:**
```
                  baseMean log2FoldChange     lfcSE      stat        pvalue          padj
ENSG00000152583   997.4398       4.313962 0.1721338  25.06168 2.378516e-138 4.341762e-134
ENSG00000165995   495.0929       3.186823 0.1281566  24.86664 3.097034e-136 2.826662e-132
ENSG00000101347 12703.3871       3.618734 0.1489434  24.29603 4.296068e-130 2.614013e-126
ENSG00000120129  3442.7126       2.871488 0.1182491  24.28339 5.823485e-130 2.657512e-126
ENSG00000189221  2341.3340       3.230395 0.1366745  23.63569 3.024222e-123 1.104082e-119
```

## 2. Medium Task ([medium.py](https://github.com/mhmohona/rosetta_gsoc_tests/blob/main/medium.py))

A pure Python script mapping the identical logic from `easy.R` into Python through `rpy2.robjects`. We extract the R dataframe natively into a classic Pandas DataFrame using `conversion.localconverter()` bounds, fully eliminating manual string bindings. 

**Output:**
```
                     baseMean  log2FoldChange     lfcSE       stat         pvalue           padj
ENSG00000152583    997.439773       -4.574919  0.184056 -24.856111  2.220933e-136  4.001677e-132
ENSG00000165995    495.092907       -3.291062  0.133174 -24.712551  7.839410e-135  7.062524e-131
ENSG00000120129   3409.029375       -2.947810  0.121438 -24.274258  3.666925e-130  2.202355e-126
ENSG00000101347  12703.387062       -3.766995  0.155438 -24.234715  9.583815e-130  4.317029e-126
ENSG00000189221   2341.767253       -3.353580  0.141782 -23.653014  1.098955e-123  3.960194e-120
```

## 3. Hard Task ([hard.py](https://github.com/mhmohona/rosetta_gsoc_tests/blob/main/medium.py))

Implementation of the programmatic `deseq2(counts, metadata, design)` wrapper pipeline.

This library strictly evaluates the incoming Pandas arrays, casting floats safely, preventing missing identifiers, and halting on illegal numeric layouts like negatives. Rather than trying to parse complex multivariant mathematical design formulas internally, we dynamically push the Pandas matrices to `rpy2.robjects` and utilize native Bioconductor `stats.as_formula()` algorithms on your `~ condition` string, guaranteeing safety around sophisticated interaction vectors without custom parser risks!

The repository also includes [test_hard.py](https://github.com/mhmohona/rosetta_gsoc_tests/blob/main/test_hard.py), implementing four rigorous internal boundary Pytest validations:
* [test_deseq2_runs_successfully](https://github.com/mhmohona/rosetta_gsoc_tests/blob/7fa13e98c5c1bb2eb331db7e456aca099f62063e/test_hard.py#L25)
* [test_input_validation_empty](https://github.com/mhmohona/rosetta_gsoc_tests/blob/7fa13e98c5c1bb2eb331db7e456aca099f62063e/test_hard.py#L36)
* [test_input_validation_matching_names](https://github.com/mhmohona/rosetta_gsoc_tests/blob/7fa13e98c5c1bb2eb331db7e456aca099f62063e/test_hard.py#L46)
* [test_input_validation_non_negative](https://github.com/mhmohona/rosetta_gsoc_tests/blob/7fa13e98c5c1bb2eb331db7e456aca099f62063e/test_hard.py#L53)

You can verify the tests utilizing Pytest:
```bash
pytest test_hard.py
```
