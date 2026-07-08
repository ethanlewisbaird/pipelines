import scanpy as sc
import pandas as pd
from scTenifold import scTenifoldKnk

# Load data
count_matrix = pd.read_csv('/data/ebaird/scRNAseq/SCENTINELsep24/composition_DEG_signatures/count_matrix.csv', index_col=0)

# # Create AnnData object
# adata = sc.AnnData(X=count_matrix.values,  # Transpose if needed
#                    var=pd.DataFrame(index=count_matrix.index),  # genes
#                    obs=pd.DataFrame(index=count_matrix.columns)) # cells

sc_tenifold = scTenifoldKnk(data=count_matrix,
                            ko_method="default",
                            ko_genes=["pros"],
                            qc_kws={"min_lib_size": 10})
result = sc_tenifold.build()

# Save results
sc_tenifold.save("./sctenifold_results")