import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

def model_a_activation_probability(stiffness, E0, k):
    return 1 / (1 + np.exp(-k * (stiffness - E0)))

def DNMT3A_active(t, decay):
    return np.exp(-decay * t)

def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

def methylation_ode(t, m, kme, kde, decay):
    dnmt3a = DNMT3A_active(t, decay)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]

def model_b_final_methylation(start_m, kme, kde, decay=0.3):
    t_span = (0, 96)
    sol = solve_ivp(methylation_ode, t_span, [start_m], args=(kme, kde, decay), method='RK45')
    return sol.y[0][-1]

GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify_mechanism(cna_status):
    cna_status = str(cna_status)
    for cat in GENETIC_LOSS_CATEGORIES:
        if cat in cna_status:
            return "genetic"
    return "epigenetic"

data_folder = "../data-acquisition/"
all_files = {
    "PAAD_CDKN2A": "cdkn2a_methylation_expression.tsv", "PAAD_SMAD4": "smad4_methylation_expression.tsv",
    "PAAD_TP53": "tp53_methylation_expression.tsv", "PAAD_BRCA1": "brca1_methylation_expression.tsv",
    "PAAD_BRCA2": "brca2_methylation_expression.tsv", "PAAD_MLH1": "mlh1_methylation_expression.tsv",
    "TNBC_CDKN2A": "tnbc_cdkn2a_methylation_expression.tsv", "TNBC_SMAD4": "tnbc_smad4_methylation_expression.tsv",
    "TNBC_TP53": "tnbc_tp53_methylation_expression.tsv", "TNBC_BRCA1": "tnbc_brca1_methylation_expression.tsv",
    "TNBC_BRCA2": "tnbc_brca2_methylation_expression.tsv", "TNBC_MLH1": "tnbc_mlh1_methylation_expression.tsv",
    "LIHC_CDKN2A": "lihc_cdkn2a_methylation_expression.tsv", "LIHC_SMAD4": "lihc_smad4_methylation_expression.tsv",
    "LIHC_TP53": "lihc_tp53_methylation_expression.tsv", "LIHC_BRCA1": "lihc_brca1_methylation_expression.tsv",
    "LIHC_BRCA2": "lihc_brca2_methylation_expression.tsv", "LIHC_MLH1": "lihc_mlh1_methylation_expression.tsv",
    "LUAD_CDKN2A": "luad_cdkn2a_methylation_expression.tsv", "LUAD_SMAD4": "luad_smad4_methylation_expression.tsv",
    "LUAD_TP53": "luad_tp53_methylation_expression.tsv", "LUAD_BRCA1": "luad_brca1_methylation_expression.tsv",
    "LUAD_BRCA2": "luad_brca2_methylation_expression.tsv", "LUAD_MLH1": "luad_mlh1_methylation_expression.tsv",
}

all_meaningful_values = []
for label, filename in all_files.items():
    path = data_folder + filename
    df = pd.read_csv(path, sep="\t")
    meth_col = [c for c in df.columns if "Methylation" in c][0]
    df_clean = df.dropna(subset=[meth_col]).copy()
    df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify_mechanism)
    epi = df_clean[df_clean["mechanism"] == "epigenetic"]
    all_meaningful_values.extend(epi[epi[meth_col] >= 0.3][meth_col].tolist())

print(f"Pooled meaningfully-silenced patients: {len(all_meaningful_values)}\n")

REACTIVATION_THRESHOLD = 0.3

def run_scenario(E0, k, kme, kde, assumed_stiffness):
    activation_prob = model_a_activation_probability(assumed_stiffness, E0, k)
    if activation_prob < 0.5:
        return 0.0
    responsive = sum(1 for m0 in all_meaningful_values if model_b_final_methylation(m0, kme, kde) < REACTIVATION_THRESHOLD)
    return (responsive / len(all_meaningful_values)) * 100

# ── RETUNED: weaker Wave 2 relative to Wave 1 (kde lowered, kme raised) ─────
# Old: kme=0.15, kde=0.4 (TET1 way overpowers DNMT3A)
# New: kme=0.25, kde=0.18 (closer fight - some patients won't fully clear by 96hrs)
print("=== Testing several kme/kde combinations to find a realistic mix ===\n")
test_combos = [
    {"kme": 0.15, "kde": 0.40, "label": "ORIGINAL (too easy)"},
    {"kme": 0.20, "kde": 0.25, "label": "Try 1"},
    {"kme": 0.25, "kde": 0.18, "label": "Try 2"},
    {"kme": 0.30, "kde": 0.12, "label": "Try 3 (harder)"},
    {"kme": 0.35, "kde": 0.08, "label": "Try 4 (much harder)"},
]

base = {"E0": 3.5, "k": 0.8, "assumed_stiffness": 4.5}
for combo in test_combos:
    result = run_scenario(E0=base["E0"], k=base["k"], kme=combo["kme"], kde=combo["kde"], assumed_stiffness=base["assumed_stiffness"])
    print(f"{combo['label']:<22} kme={combo['kme']}, kde={combo['kde']}  ->  {result:.1f}% responsive")

print("\n=== Full sensitivity analysis using the best realistic combo (Try 3) ===\n")
baseline = {"E0": 3.5, "k": 0.8, "kme": 0.30, "kde": 0.12, "assumed_stiffness": 4.5}
baseline_result = run_scenario(**baseline)
print(f"BASELINE % responsive: {baseline_result:.1f}%\n")

print(f"{'Parameter':<10} {'Change':<10} {'New value':<12} {'% responsive':<14} {'Shift'}")
print("-" * 70)
for param in ["E0", "k", "kme", "kde"]:
    for pct in [-0.3, -0.2, 0.2, 0.3]:
        test_params = baseline.copy()
        test_params[param] = baseline[param] * (1 + pct)
        result = run_scenario(**test_params)
        shift = result - baseline_result
        print(f"{param:<10} {pct*100:+.0f}%{'':<6} {test_params[param]:<12.3f} {result:<14.1f} {shift:+.1f} pts")
