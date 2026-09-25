import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

# ── Model A: Mechanosensing Trigger ──────────────────────────────────────────
# Sigmoid converts tissue stiffness (kPa) into a release probability (0→1)
# UPDATED: E0/k now based on real MR elastography literature values,
# replacing earlier placeholder estimates. LUAD reuses LIHC's calibration
# since no solid MRE tumor-vs-normal pair was found in the literature for lung
# (stated limitation, not guessed silently).

def activation_probability(E, E0, k):
    """
    E  = tissue stiffness in kPa (the input)
    E0 = threshold stiffness where release is 50% likely
    k  = steepness of the transition (higher = sharper switch)
    """
    return 1 / (1 + np.exp(-k * (E - E0)))

# ── Cancer-specific parameters (REAL literature values, MRE-based) ──────────
# Healthy/tumor ranges and E0 derived from MR elastography studies.
# k chosen to reflect how WIDE the healthy/tumor gap is per cancer:
#   - Pancreatic: gap is small (~1 to ~3-6 kPa) -> shallow k (gradual, overlapping)
#   - Liver: moderate gap (~2-3 to ~7-10 kPa), confounded by cirrhosis -> moderate k
#   - Breast: large gap (~2-2.5 to ~15.9 kPa) -> steep k (clean separation)
#   - Lung: NO independent MRE data found -> reuses LIHC's E0/k as a stated,
#     documented assumption rather than an undocumented guess.

cancers = {
    "Pancreatic (PAAD)": {
        "E0": 3.5, "k": 0.8, "color": "#2a78d6",
        "healthy_range": (0.7, 1.3),   # ~1 kPa normal pancreas (MRE)
        "tumor_range": (3, 6),         # ~3-6 kPa PDAC tumor (MRE)
        "note": "Weak separation - healthy/tumor stiffness nearly overlap"
    },
    "Hepatocellular (LIHC)": {
        "E0": 5.5, "k": 0.6, "color": "#eda100",
        "healthy_range": (2, 3),       # cirrhotic background liver (MRE)
        "tumor_range": (7, 10),        # HCC tumor (MRE)
        "note": "Moderate separation - confounded by background cirrhosis"
    },
    "Triple-neg Breast (TNBC)": {
        "E0": 9.0, "k": 0.35, "color": "#1baf7a",
        "healthy_range": (2, 2.5),     # normal breast parenchyma (MRE)
        "tumor_range": (12, 20),       # malignant breast tumor (MRE, ~15.9 kPa median)
        "note": "Strong separation - clean discrimination"
    },
    "Lung Adeno (LUAD)": {
        "E0": 5.5, "k": 0.6, "color": "#333333",
        "healthy_range": (2, 3),
        "tumor_range": (7, 10),
        "note": "Reused LIHC calibration - no MRE tumor stiffness data found for lung (stated limitation)"
    },
}

# ── Plot ──────────────────────────────────────────────────────────────────────
E_values = np.linspace(0, 25, 500)   # narrower range - real values are much lower than old 0-60 placeholder

fig, ax = plt.subplots(figsize=(10, 6))

for name, params in cancers.items():
    P = activation_probability(E_values, params["E0"], params["k"])
    is_luad = "LUAD" in name
    linestyle = "--" if is_luad else "-"
    lw = 3.5 if is_luad else 2.5
    zorder = 10 if is_luad else 1
    ax.plot(E_values, P, label=name, color=params["color"], linewidth=lw, linestyle=linestyle, zorder=zorder)
    if not is_luad:
        ax.axvspan(params["healthy_range"][0], params["healthy_range"][1],
                   alpha=0.10, color=params["color"])
        ax.axvspan(params["tumor_range"][0], params["tumor_range"][1],
                   alpha=0.06, color=params["color"])

ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax.text(25.3, 0.5, "50%", va="center", fontsize=10, color="gray")

ax.set_xlabel("Tissue stiffness (kPa) — MR Elastography scale", fontsize=13)
ax.set_ylabel("MENR activation probability", fontsize=13)
ax.set_title("Model A — Mechanosensing Trigger (corrected with real MRE literature values)\n"
             "Note: LUAD (dashed, black) reuses LIHC calibration - no solid MRE tumor stiffness data found in literature",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="center right")
ax.set_xlim(0, 25)
ax.set_ylim(-0.02, 1.05)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../figures/model_a_output_v2.png", dpi=150, bbox_inches='tight')
print("Saved! Check figures/model_a_output_v2.png")

# Print a summary table so the differences are obvious in terminal too
print("\n--- Calibration summary (real MRE literature values) ---")
for name, p in cancers.items():
    print(f"{name}: E0={p['E0']} kPa, k={p['k']}, healthy={p['healthy_range']}, tumor={p['tumor_range']}")
    print(f"  -> {p['note']}")
