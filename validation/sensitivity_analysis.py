import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

# ── Same core functions as the digital twin, but parameterized for testing ──
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

def model_b_final_methylation(start_m, kme, kde, decay=0.3, threshold=0.3):
    t_span = (0, 96)
    sol = solve_ivp(methylation_ode, t_span, [start_m], args=(kme, kde, decay), method='RK45')
    final_m = sol.y[0][-1]
    return final_m < threshold

GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify_mechanism(cna_status):
    cna_status = str(cna_status)
    for cat in GENETIC_LOSS_CATEGORIES:
        if cat in cna_status:
            return "genetic"
    return "epigenetic"

# ── Load PANCREATIC CDKN2A data as our test case (smaller, faster to test) ──
data_path = "../data-acquisition/cdkn2a_methylation_expression.tsv"
df = pd.read_csv(data_path, sep="\t")
meth_col = [c for c in df.columns if "Methylation" in c][0]
df_clean = df.dropna(subset=[meth_col]).copy()
df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify_mechanism)
epi_patients = df_clean[df_clean["mechanism"] == "epigenetic"]

# ── Baseline parameters (what we've been using) ─────────────────────────────
baseline = {"E0": 3.5, "k": 0.8, "kme": 0.15, "kde": 0.4, "assumed_stiffness": 4.5}

def run_scenario(E0, k, kme, kde, assumed_stiffness):
    activation_prob = model_a_activation_probability(assumed_stiffness, E0, k)
    if activation_prob < 0.5:
        return 0.0  # MENR doesn't activate, nobody responds
    responsive_count = 0
    for _, row in epi_patients.iterrows():
        if model_b_final_methylation(row[meth_col], kme, kde):
            responsive_count += 1
    return (responsive_count / len(epi_patients)) * 100

# ── Run baseline ──────────────────────────────────────────────────────────
baseline_result = run_scenario(**baseline)
print(f"BASELINE % responsive (epigenetic patients only): {baseline_result:.1f}%\n")

# ── Vary each parameter by +/-20% and +/-30%, one at a time ────────────────
print(f"{'Parameter':<10} {'Change':<10} {'New value':<12} {'% responsive':<14} {'Shift from baseline'}")
print("-" * 70)

params_to_test = ["E0", "k", "kme", "kde"]
percent_changes = [-0.3, -0.2, 0.2, 0.3]

for param in params_to_test:
    for pct in percent_changes:
        test_params = baseline.copy()
        new_value = baseline[param] * (1 + pct)
        test_params[param] = new_value
        result = run_scenario(**test_params)
        shift = result - baseline_result
        print(f"{param:<10} {pct*100:+.0f}%{'':<6} {new_value:<12.3f} {result:<14.1f} {shift:+.1f} pts")
