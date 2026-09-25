import pandas as pd

genes = {
    "CDKN2A": "tnbc_cdkn2a_methylation_expression.tsv",
    "SMAD4": "tnbc_smad4_methylation_expression.tsv",
    "TP53": "tnbc_tp53_methylation_expression.tsv",
    "BRCA1": "tnbc_brca1_methylation_expression.tsv",
    "BRCA2": "tnbc_brca2_methylation_expression.tsv",
    "MLH1": "tnbc_mlh1_methylation_expression.tsv",
}

data_folder = "../data-acquisition/"
GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]

def classify(cna_status):
    cna_status = str(cna_status)
    for category in GENETIC_LOSS_CATEGORIES:
        if category in cna_status:
            return "genetic"
    return "epigenetic"

print("BREAST CANCER (full cohort, not yet filtered to TNBC) - gene breakdown\n")
print(f"{'Gene':<10} {'Total':<8} {'Genetic':<10} {'Epigenetic':<12} {'% treatable'}")
print("-" * 60)

for gene, filename in genes.items():
    path = data_folder + filename
    df = pd.read_csv(path, sep="\t")
    meth_col = [c for c in df.columns if "Methylation" in c][0]
    df_clean = df.dropna(subset=[meth_col]).copy()
    df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify)

    total = len(df_clean)
    genetic_count = (df_clean["mechanism"] == "genetic").sum()
    epigenetic_count = (df_clean["mechanism"] == "epigenetic").sum()
    pct_treatable = (epigenetic_count / total * 100) if total > 0 else 0

    print(f"{gene:<10} {total:<8} {genetic_count:<10} {epigenetic_count:<12} {pct_treatable:.1f}%")
