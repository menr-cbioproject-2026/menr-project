import pandas as pd

# All 6 genes we pulled data for, and their file names
genes = {
    "CDKN2A": "cdkn2a_methylation_expression.tsv",
    "SMAD4": "smad4_methylation_expression.tsv",
    "TP53": "tp53_methylation_expression.tsv",
    "BRCA1": "brca1_methylation_expression.tsv",
    "BRCA2": "brca2_methylation_expression.tsv",
    "MLH1": "mlh1_methylation_expression.tsv",
}

data_folder = "../data-acquisition/"

# Same rule as before: deletion = gene is gone = MENR can't fix it
GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]

def classify(cna_status):
    cna_status = str(cna_status)
    for category in GENETIC_LOSS_CATEGORIES:
        if category in cna_status:
            return "genetic"
    return "epigenetic"

print("Gene-by-gene breakdown: how many patients MENR could even theoretically help\n")
print(f"{'Gene':<10} {'Total':<8} {'Genetic (cant help)':<22} {'Epigenetic (could help)':<24} {'% treatable'}")
print("-" * 85)

results = []

for gene, filename in genes.items():
    path = data_folder + filename
    df = pd.read_csv(path, sep="\t")

    meth_col = [c for c in df.columns if "Methylation" in c][0]
    df_clean = df.dropna(subset=[meth_col])

    df_clean = df_clean.copy()
    df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify)

    total = len(df_clean)
    genetic_count = (df_clean["mechanism"] == "genetic").sum()
    epigenetic_count = (df_clean["mechanism"] == "epigenetic").sum()
    pct_treatable = (epigenetic_count / total * 100) if total > 0 else 0

    results.append({
        "gene": gene, "total": total,
        "genetic": genetic_count, "epigenetic": epigenetic_count,
        "pct_treatable": pct_treatable
    })

    print(f"{gene:<10} {total:<8} {genetic_count:<22} {epigenetic_count:<24} {pct_treatable:.1f}%")

print("\n--- In plain terms ---")
for r in results:
    print(f"{r['gene']}: out of {r['total']} patients, only {r['pct_treatable']:.0f}% even have the gene silenced in a way MENR could fix.")
