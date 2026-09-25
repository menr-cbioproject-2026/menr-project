import pandas as pd

data_path = "../data-acquisition/cdkn2a_methylation_expression.tsv"
df = pd.read_csv(data_path, sep="\t")
meth_col = [c for c in df.columns if "Methylation" in c][0]
df_clean = df.dropna(subset=[meth_col])

GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify(cna_status):
    cna_status = str(cna_status)
    for cat in GENETIC_LOSS_CATEGORIES:
        if cat in cna_status:
            return "genetic"
    return "epigenetic"

df_clean = df_clean.copy()
df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify)
epi = df_clean[df_clean["mechanism"] == "epigenetic"]

print("Starting methylation values - epigenetic patients only:")
print(epi[meth_col].describe())
print(f"\nLowest 5 values:\n{epi[meth_col].nsmallest(5).values}")
print(f"\nHighest 5 values:\n{epi[meth_col].nlargest(5).values}")
