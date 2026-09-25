import pandas as pd

clin = pd.read_csv('data-acquisition/lihc_tcga_pan_can_atlas_2018_clinical_data.tsv', sep='\t')
meth = pd.read_csv('data-acquisition/lihc_rassf1_methylation_expression.tsv', sep='\t')

stage_col = 'American Joint Committee on Cancer Tumor Stage Code'
print("Stage col exists:", stage_col in clin.columns)

clin_sub = clin[['Sample ID', stage_col]].rename(columns={stage_col: 'stage'})
meth_col = [c for c in meth.columns if 'Methylation' in c][0]
meth_sub = meth[['Sample Id', meth_col]].rename(columns={'Sample Id': 'Sample ID', meth_col: 'methylation'})

print("clin_sub shape:", clin_sub.shape)
print("meth_sub shape:", meth_sub.shape)
print("clin_sub Sample ID dtype:", clin_sub['Sample ID'].dtype)
print("meth_sub Sample ID dtype:", meth_sub['Sample ID'].dtype)

merged = pd.merge(meth_sub, clin_sub, on='Sample ID')
print("MERGED shape:", merged.shape)
print(merged.head())
