import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify(cna_status):
    cna_status = str(cna_status)
    for cat in GENETIC_LOSS_CATEGORIES:
        if cat in cna_status:
            return "genetic"
    return "epigenetic"

def DNMT3A_active(t, natural_level=0.05, inhibitor_strength=0.85, inhibitor_decay=0.12):
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)

def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

def methylation_ode(t, m, kme=0.20, kde=0.30):
    return [kme * (1 - m) * DNMT3A_active(t) - kde * m * TET1_active(t)]

def model_b_outcome(start_m):
    t_eval = np.linspace(0, 96, 100)
    sol = solve_ivp(methylation_ode, (0, 96), [start_m], t_eval=t_eval, method='RK45')
    final_m = sol.y[0][-1]
    pct_drop = (start_m - final_m) / start_m * 100 if start_m > 0 else 0
    return pct_drop >= 40

data_folder = "../data-acquisition/"

# 5 genes x 4 cancers - mix of strong (CDKN2A, RASSF1, APC) and weak (SMAD4, TP53)
# original panel genes, shown side by side for an honest full picture
gene_files = {
    "PAAD": {
        "CDKN2A": "cdkn2a_methylation_expression.tsv",
        "RASSF1": "paad_rassf1_methylation_expression.tsv",
        "APC": "paad_apc_methylation_expression.tsv",
        "SMAD4": "smad4_methylation_expression.tsv",
        "TP53": "tp53_methylation_expression.tsv",
    },
    "TNBC": {
        "CDKN2A": "tnbc_cdkn2a_methylation_expression.tsv",
        "RASSF1": "tnbc_rassf1_methylation_expression.tsv",
        "APC": "tnbc_apc_methylation_expression.tsv",
        "SMAD4": "tnbc_smad4_methylation_expression.tsv",
        "TP53": "tnbc_tp53_methylation_expression.tsv",
    },
    "LIHC": {
        "CDKN2A": "lihc_cdkn2a_methylation_expression.tsv",
        "RASSF1": "lihc_rassf1_methylation_expression.tsv",
        "APC": "lihc_apc_methylation_expression.tsv",
        "SMAD4": "lihc_smad4_methylation_expression.tsv",
        "TP53": "lihc_tp53_methylation_expression.tsv",
    },
    "LUAD": {
        "CDKN2A": "luad_cdkn2a_methylation_expression.tsv",
        "RASSF1": "luad_rassf1_methylation_expression.tsv",
        "APC": "luad_apc_methylation_expression.tsv",
        "SMAD4": "luad_smad4_methylation_expression.tsv",
        "TP53": "luad_tp53_methylation_expression.tsv",
    },
}

rows = []
for cancer, genes in gene_files.items():
    for gene, filename in genes.items():
        df = pd.read_csv(data_folder + filename, sep="\t")
        meth_col = [c for c in df.columns if "Methylation" in c][0]
        df_clean = df.dropna(subset=[meth_col]).copy()
        df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify)

        total = len(df_clean)
        epi = df_clean[df_clean["mechanism"] == "epigenetic"]
        heavily_silenced = epi[epi[meth_col] >= 0.3]

        responsive_count = sum(1 for _, row in heavily_silenced.iterrows() if model_b_outcome(row[meth_col]))
        n_eligible = len(heavily_silenced)
        pct_resp = (responsive_count / n_eligible * 100) if n_eligible > 0 else 0
        pct_eligible = (n_eligible / total * 100) if total > 0 else 0

        rows.append({
            "cancer": cancer, "gene": gene, "total_patients": total,
            "pct_eligible_heavily_silenced": round(pct_eligible, 1),
            "n_eligible": n_eligible,
            "pct_predicted_responsive": round(pct_resp, 1),
        })

results = pd.DataFrame(rows)
results.to_csv("../digital-twin/master_comparison_table.csv", index=False)

print(results.to_string(index=False))
print("\nSaved to master_comparison_table.csv")

print("\n--- Summary: average % eligible by gene (across all 4 cancers) ---")
print(results.groupby("gene")["pct_eligible_heavily_silenced"].mean().round(1))

print("\n--- Summary: average % predicted responsive by gene (among eligible) ---")
print(results.groupby("gene")["pct_predicted_responsive"].mean().round(1))
