import pandas as pd

data_folder = "../data-acquisition/"

all_files = {
    "PAAD_CDKN2A": "cdkn2a_methylation_expression.tsv", "PAAD_SMAD4": "smad4_methylation_expression.tsv",
    "PAAD_TP53": "tp53_methylation_expression.tsv", "PAAD_BRCA1": "brca1_methylation_expression.tsv",
    "PAAD_BRCA2": "brca2_methylation_expression.tsv", "PAAD_MLH1": "mlh1_methylation_expression.tsv",
    "TNBC_CDKN2A": "tnbc_cdkn2a_methylation_expression.tsv", "TNBC_SMAD4": "tnbc_smad4_methylation_expression.tsv",
    "TNBC_TP53": "tnbc_tp53_methylation_expression.tsv", "TNBC_BRCA1": "tnbc_brca1_methylation_expression.tsv",
    "TNBC_BRCA2": "tnbc_brca2_methylation_expression.tsv", "TNBC_MLH1": "tnbc_mlh1_methylation_expression.tsv",
}

print(f"{'File':<16}{'Count':<8}{'Mean':<8}{'Max':<8}{'% above 0.55 (can reach 40% recovery)'}")
for label, filename in all_files.items():
    df = pd.read_csv(data_folder + filename, sep="\t")
    meth_col = [c for c in df.columns if "Methylation" in c][0]
    vals = df[meth_col].dropna()
    pct_above = (vals >= 0.55).mean() * 100
    print(f"{label:<16}{len(vals):<8}{vals.mean():<8.3f}{vals.max():<8.3f}{pct_above:.1f}%")
