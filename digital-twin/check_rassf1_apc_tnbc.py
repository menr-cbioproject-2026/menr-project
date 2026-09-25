import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

genes = {
    "RASSF1": "tnbc_rassf1_methylation_expression.tsv",
    "APC": "tnbc_apc_methylation_expression.tsv",
}

data_folder = "../data-acquisition/"
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
    return final_m, pct_drop, pct_drop >= 40

print("BREAST (TNBC, full cohort) - RASSF1 and APC\n")
for gene, filename in genes.items():
    df = pd.read_csv(data_folder + filename, sep="\t")
    meth_col = [c for c in df.columns if "Methylation" in c][0]
    df_clean = df.dropna(subset=[meth_col]).copy()
    df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify)

    total = len(df_clean)
    vals = df_clean[meth_col]
    pct_above_55 = (vals >= 0.55).mean() * 100

    epi = df_clean[df_clean["mechanism"] == "epigenetic"]
    heavily_silenced = epi[epi[meth_col] >= 0.3]

    responsive_count = 0
    for _, row in heavily_silenced.iterrows():
        _, _, reactivated = model_b_outcome(row[meth_col])
        if reactivated:
            responsive_count += 1

    print(f"--- {gene} ---")
    print(f"Total patients: {total}")
    print(f"Mean methylation: {vals.mean():.3f}, Max: {vals.max():.3f}")
    print(f"% of patients above 0.55 methylation: {pct_above_55:.1f}%")
    print(f"Epigenetic (not deleted) patients: {len(epi)}")
    print(f"Heavily silenced (>=0.3) AND epigenetic: {len(heavily_silenced)}")
    pct_resp = responsive_count/len(heavily_silenced)*100 if len(heavily_silenced)>0 else 0
    print(f"Of those, predicted responsive to MENR: {responsive_count} ({pct_resp:.1f}%)")
    print()
