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

data_path = "../data-acquisition/cdkn2a_methylation_expression.tsv"
df = pd.read_csv(data_path, sep="\t")
meth_col = [c for c in df.columns if "Methylation" in c][0]
df_clean = df.dropna(subset=[meth_col]).copy()
df_clean["mechanism"] = df_clean["Copy Number Alterations"].apply(classify_mechanism)
epi_patients = df_clean[df_clean["mechanism"] == "epigenetic"]

# ── NEW: only patients who START genuinely silenced (above 0.3) count ──────
# Below 0.3 = gene was never really "switched off" - nothing for MENR to fix
MEANINGFULLY_SILENCED_THRESHOLD = 0.3
testable_patients = epi_patients[epi_patients[meth_col] >= MEANINGFULLY_SILENCED_THRESHOLD]
print(f"Total epigenetic patients: {len(epi_patients)}")
print(f"Patients actually meaningfully silenced (start >= 0.3): {len(testable_patients)}")
print(f"These are the only ones where 'success' is a meaningful claim.\n")

baseline = {"E0": 3.5, "k": 0.8, "kme": 0.15, "kde": 0.4, "assumed_stiffness": 4.5}
REACTIVATION_THRESHOLD = 0.3

def run_scenario(E0, k, kme, kde, assumed_stiffness):
    activation_prob = model_a_activation_probability(assumed_stiffness, E0, k)
    if activation_prob < 0.5 or len(testable_patients) == 0:
        return 0.0
    responsive_count = 0
    for _, row in testable_patients.iterrows():
        final_m = model_b_final_methylation(row[meth_col], kme, kde)
        if final_m < REACTIVATION_THRESHOLD:
            responsive_count += 1
    return (responsive_count / len(testable_patients)) * 100

baseline_result = run_scenario(**baseline)
print(f"BASELINE % responsive (meaningfully-silenced patients only): {baseline_result:.1f}%\n")

print(f"{'Parameter':<10} {'Change':<10} {'New value':<12} {'% responsive':<14} {'Shift'}")
print("-" * 70)
for param in ["E0", "k", "kme", "kde"]:
    for pct in [-0.3, -0.2, 0.2, 0.3]:
        test_params = baseline.copy()
        test_params[param] = baseline[param] * (1 + pct)
        result = run_scenario(**test_params)
        shift = result - baseline_result
        print(f"{param:<10} {pct*100:+.0f}%{'':<6} {test_params[param]:<12.3f} {result:<14.1f} {shift:+.1f} pts")
