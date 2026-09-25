import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

# ── PART 1: Model A - stiffness activation, real literature values per cancer ─
cancer_stiffness_params = {
    "PAAD":  {"E0": 3.5, "k": 0.8},
    "LIHC":  {"E0": 5.5, "k": 0.6},
    "TNBC":  {"E0": 9.0, "k": 0.35},
    "LUAD":  {"E0": 5.5, "k": 0.6},  # reusing liver's k - no solid lung MRE pair found (stated limitation)
}

def model_a_activation_probability(cancer_type, assumed_stiffness_kpa):
    params = cancer_stiffness_params[cancer_type]
    E0, k = params["E0"], params["k"]
    return 1 / (1 + np.exp(-k * (assumed_stiffness_kpa - E0)))

assumed_tumor_stiffness = {
    "PAAD": 4.5, "LIHC": 8.5, "TNBC": 16.0, "LUAD": 8.5,
}

# ── PART 2: Genetic vs epigenetic classification ───────────────────────────
GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify_mechanism(cna_status):
    cna_status = str(cna_status)
    for category in GENETIC_LOSS_CATEGORIES:
        if category in cna_status:
            return "genetic"
    return "epigenetic"

# ── PART 3: CORRECTED Model B kinetics ──────────────────────────────────────
# DNMT3A: ongoing background presence (Section 2.3), Wave 1 inhibitor's
# effect on it fades over time - DNMT3A itself does not die out.
# NOTE: natural_level, inhibitor_strength/decay, kme, kde below are TUNED
# FIT PARAMETERS (no published in-vitro rate exists for this specific
# reaction) - chosen to produce a biologically plausible mixed responder/
# non-responder pattern, per Master Doc Section 4.2.1's allowance for
# stated-uncertainty fit parameters. This must be disclosed as such in
# the report (Section 13.3).
def DNMT3A_active(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)

def TET1_active(t, delay=12, rise=0.3, clearance=0.08):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

def methylation_ode(t, m, kme=0.15, kde=0.40):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]

RECOVERY_FRACTION = 0.40  # relative recovery threshold (Section 4.2.2 - "expression-recovery threshold")

def model_b_outcome(start_m):
    t_span = (0, 48)
    t_eval = np.linspace(0, 48, 100)
    solution = solve_ivp(methylation_ode, t_span, [start_m], t_eval=t_eval, method='RK45')
    final_m = solution.y[0][-1]
    pct_drop = (start_m - final_m) / start_m * 100 if start_m > 0 else 0
    reactivated = pct_drop >= (RECOVERY_FRACTION * 100)
    return final_m, pct_drop, reactivated

# ── PART 4: Run across all 4 cancers x 6 genes x real patients ─────────────
data_folder = "../data-acquisition/"

cancer_gene_files = {
    "PAAD": {
        "CDKN2A": "cdkn2a_methylation_expression.tsv", "SMAD4": "smad4_methylation_expression.tsv",
        "TP53": "tp53_methylation_expression.tsv", "BRCA1": "brca1_methylation_expression.tsv",
        "BRCA2": "brca2_methylation_expression.tsv", "MLH1": "mlh1_methylation_expression.tsv", "RASSF1": "paad_rassf1_methylation_expression.tsv", "APC": "paad_apc_methylation_expression.tsv",
    },
    "TNBC": {
        "CDKN2A": "tnbc_cdkn2a_methylation_expression.tsv", "SMAD4": "tnbc_smad4_methylation_expression.tsv",
        "TP53": "tnbc_tp53_methylation_expression.tsv", "BRCA1": "tnbc_brca1_methylation_expression.tsv",
        "BRCA2": "tnbc_brca2_methylation_expression.tsv", "MLH1": "tnbc_mlh1_methylation_expression.tsv", "RASSF1": "tnbc_rassf1_methylation_expression.tsv", "APC": "tnbc_apc_methylation_expression.tsv",
    },
    "LIHC": {
        "CDKN2A": "lihc_cdkn2a_methylation_expression.tsv", "SMAD4": "lihc_smad4_methylation_expression.tsv",
        "TP53": "lihc_tp53_methylation_expression.tsv", "BRCA1": "lihc_brca1_methylation_expression.tsv",
        "BRCA2": "lihc_brca2_methylation_expression.tsv", "MLH1": "lihc_mlh1_methylation_expression.tsv", "RASSF1": "lihc_rassf1_methylation_expression.tsv", "APC": "lihc_apc_methylation_expression.tsv",
    },
    "LUAD": {
        "CDKN2A": "luad_cdkn2a_methylation_expression.tsv", "SMAD4": "luad_smad4_methylation_expression.tsv",
        "TP53": "luad_tp53_methylation_expression.tsv", "BRCA1": "luad_brca1_methylation_expression.tsv",
        "BRCA2": "luad_brca2_methylation_expression.tsv", "MLH1": "luad_mlh1_methylation_expression.tsv", "RASSF1": "luad_rassf1_methylation_expression.tsv", "APC": "luad_apc_methylation_expression.tsv",
    },
}

all_results = []

for cancer_type, gene_files in cancer_gene_files.items():
    activation_prob = model_a_activation_probability(cancer_type, assumed_tumor_stiffness[cancer_type])

    for gene, filename in gene_files.items():
        path = data_folder + filename
        df = pd.read_csv(path, sep="\t")
        meth_col = [c for c in df.columns if "Methylation" in c][0]
        df_clean = df.dropna(subset=[meth_col]).copy()

        for _, row in df_clean.iterrows():
            patient_id = row["Sample Id"]
            cna_status = row["Copy Number Alterations"]
            start_m = row[meth_col]
            mechanism = classify_mechanism(cna_status)

            if mechanism == "genetic":
                predicted_responsive = False
                final_m, pct_drop = None, None
                reason = "Genetic deletion - epigenetic payload cannot apply"
            elif activation_prob < 0.5:
                predicted_responsive = False
                final_m, pct_drop = None, None
                reason = f"MENR unlikely to activate (stiffness probability={activation_prob:.2f})"
            elif start_m <= 0:
                predicted_responsive = False
                final_m, pct_drop = None, None
                reason = "Starting methylation is zero - no silencing to reverse"
            else:
                final_m, pct_drop, reactivated = model_b_outcome(start_m)
                predicted_responsive = reactivated
                reason = f"Full pipeline: {pct_drop:.1f}% methylation reduction over 48hrs"

            all_results.append({
                "patient_id": patient_id,
                "cancer_type": cancer_type,
                "gene": gene,
                "cna_status": cna_status,
                "mechanism": mechanism,
                "activation_probability": round(activation_prob, 3),
                "starting_methylation": start_m,
                "final_methylation": round(final_m, 4) if final_m is not None else None,
                "pct_methylation_drop": round(pct_drop, 1) if pct_drop is not None else None,
                "predicted_responsive": predicted_responsive,
                "reason": reason,
            })

results_df = pd.DataFrame(all_results)
results_df.to_csv("../digital-twin/digital_twin_results_v2.csv", index=False)

print(f"Total patient-gene predictions made: {len(results_df)}")
print(f"\nOverall responsive rate: {results_df['predicted_responsive'].mean()*100:.1f}%")
print("\nResponsive rate by cancer type:")
print(results_df.groupby("cancer_type")["predicted_responsive"].mean() * 100)
print("\nResponsive rate by gene:")
print(results_df.groupby("gene")["predicted_responsive"].mean() * 100)
print("\nResponsive rate by silencing mechanism:")
print(results_df.groupby("mechanism")["predicted_responsive"].mean() * 100)
print("\nSaved full results to digital_twin_results_v2.csv")
