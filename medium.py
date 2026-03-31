import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def run_medium_test():
    """
    Run DESeq2 on airway and convert results to a pandas DataFrame.
    """
    ro.r('.libPaths(c("~/R/libs", .libPaths()))')
    
    importr('DESeq2')
    importr('airway')

    ro.r('''
    data(airway)
    dds <- DESeqDataSet(airway, design = ~ cell + dex)
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    df_res <- as.data.frame(resOrdered)
    ''')

    with localconverter(ro.default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(ro.r['df_res'])

    print(res_df.head(5))
    return res_df

if __name__ == "__main__":
    run_medium_test()
