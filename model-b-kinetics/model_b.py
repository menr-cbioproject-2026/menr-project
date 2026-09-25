import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ── Load REAL patient methylation data (CDKN2A, pancreatic cancer, TCGA) ────
data_path = "../data-acquisition/cdkn2a_methylation_expression.tsv"
df = pd.read_csv(data_path, sep="\t")
meth_col = [c for c in df.columns if "Methylation" in c][0]
df_clean = df.dropna(subset=[meth_col])

GENETIC_LOSS_CATEGORIES = ["Deep Deletion", "Shallow Deletion"]
def classify_silencing_mechanism(cna_status):
    cna_status = str(cna_status)
    for category in GENETIC_LOSS_CATEGORIES:
        if category in cna_status:
            return "genetic"
    return "epigenetic"

df_clean = df_clean.copy()
df_clean["silencing_mechanism"] = df_clean["Copy Number Alterations"].apply(classify_silencing_mechanism)

# ── CORRECTED Wave 1: DNMT3A keeps an ongoing background presence ──────────
# (matches Master Doc Section 2.3: "Wave 2 alone is continually counteracted
# by ongoing DNMT3A activity" - DNMT3A does NOT fully die out, only the
# Wave 1 inhibitor's effect on it fades over time)
def DNMT3A_active(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
    # natural_level raised to 0.12 (from 0.05) so background DNMT3A provides
    # meaningful competition against TET1 - tuned fit parameter, stated uncertainty
    # inhibitor_decay raised to 0.4 so Wave 1 suppression fades realistically
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)

# ── Wave 2: TET1 production over time (unchanged) ──────────────────────────
def TET1_active(t, delay=12, rise=0.3, clearance=0.08):
    # clearance raised to 0.08 (from 0.05) so TET1 fades within 48hr window
    # this creates the biologically realistic gradient: high starters succeed
    # before TET1 clears, low starters do not - tuned fit parameter
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

# ── The ODE: methylation change over time ───────────────────────────────────
def methylation_ode(t, m, kme=0.15, kde=0.40):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]

# ── CORRECTED success rule: relative recovery, not an absolute cutoff ──────
# (matches Section 4.2.2's "expression-recovery threshold" - recovery is
# inherently relative to where a patient started, not one universal number)
RECOVERY_FRACTION = 0.40  # must drop by >=40% from starting methylation to count as responsive

def simulate_patient(real_m0, mechanism):
    if mechanism == "genetic":
        return {
            "t": None, "m": None, "final_m": None, "pct_drop": None,
            "predicted_reactivated": False,
            "reason": "Gene silenced by genetic deletion - epigenetic payload cannot apply"
        }
    t_span = (0, 48)
    t_eval = np.linspace(0, 48, 200)
    solution = solve_ivp(methylation_ode, t_span, [real_m0], t_eval=t_eval, method='RK45')
    t = solution.t
    m = solution.y[0]
    final_m = m[-1]
    pct_drop = (real_m0 - final_m) / real_m0 * 100
    reactivated = pct_drop >= (RECOVERY_FRACTION * 100)
    return {
        "t": t, "m": m, "final_m": final_m, "pct_drop": pct_drop,
        "predicted_reactivated": reactivated,
        "reason": f"Epigenetic silencing - {pct_drop:.1f}% methylation reduction over 48hrs"
    }

# ── Run on first genetic patient and first epigenetic patient, for comparison ─
genetic_patients = df_clean[df_clean["silencing_mechanism"] == "genetic"]
epi_patients = df_clean[df_clean["silencing_mechanism"] == "epigenetic"]

gen_row = genetic_patients.iloc[0]
gen_result = simulate_patient(gen_row[meth_col], "genetic")
print(f"--- Genetic deletion patient: {gen_row['Sample Id']} ---")
print(f"CNA status: {gen_row['Copy Number Alterations']}")
print(f"Predicted reactivated: {gen_result['predicted_reactivated']}")
print(f"Reason: {gen_result['reason']}\n")

# Find a LOW-starting and a HIGH-starting epigenetic patient to show the new mixed behavior
epi_sorted = epi_patients.sort_values(meth_col)
low_row = epi_sorted.iloc[0]
high_row = epi_sorted.iloc[-1]

low_result = simulate_patient(low_row[meth_col], "epigenetic")
high_result = simulate_patient(high_row[meth_col], "epigenetic")

print(f"--- Low-starting epigenetic patient: {low_row['Sample Id']} ---")
print(f"Starting methylation: {low_row[meth_col]:.4f}")
print(f"Final methylation: {low_result['final_m']:.4f} ({low_result['pct_drop']:.1f}% drop)")
print(f"Predicted reactivated: {low_result['predicted_reactivated']}\n")

print(f"--- High-starting epigenetic patient: {high_row['Sample Id']} ---")
print(f"Starting methylation: {high_row[meth_col]:.4f}")
print(f"Final methylation: {high_result['final_m']:.4f} ({high_result['pct_drop']:.1f}% drop)")
print(f"Predicted reactivated: {high_result['predicted_reactivated']}")

# ── Plot: genetic case vs low-starter vs high-starter, side by side ────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

axes[0].axhline(gen_row[meth_col], color="#999999", linewidth=3)
axes[0].set_title(f"GENETIC DELETION\n{gen_row['Sample Id']}\nReactivated: {gen_result['predicted_reactivated']}", fontsize=11, color="#c0392b")
axes[0].set_ylim(0, 1); axes[0].set_xlim(0, 48)
axes[0].set_xlabel("Time (hours)"); axes[0].set_ylabel("Methylation fraction")
axes[0].grid(True, alpha=0.3)

axes[1].plot(low_result["t"], low_result["m"], color="#e67e22", linewidth=3)
axes[1].set_title(f"LOW STARTER (mild silencing)\nStart={low_row[meth_col]:.2f} -> {low_result['final_m']:.2f}\nReactivated: {low_result['predicted_reactivated']}", fontsize=11, color="#b9770e")
axes[1].set_ylim(0, 1); axes[1].set_xlim(0, 48)
axes[1].set_xlabel("Time (hours)")
axes[1].grid(True, alpha=0.3)

axes[2].plot(high_result["t"], high_result["m"], color="#1baf7a", linewidth=3)
axes[2].set_title(f"HIGH STARTER (heavy silencing)\nStart={high_row[meth_col]:.2f} -> {high_result['final_m']:.2f}\nReactivated: {high_result['predicted_reactivated']}", fontsize=11, color="#1e8449")
axes[2].set_ylim(0, 1); axes[2].set_xlim(0, 48)
axes[2].set_xlabel("Time (hours)")
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../figures/model_b_corrected.png", dpi=150, bbox_inches='tight')
print("\nSaved! Check figures/model_b_corrected.png")
