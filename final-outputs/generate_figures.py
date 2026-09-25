import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import solve_ivp

df = pd.read_csv('../digital-twin/digital_twin_results_v2.csv')
outdir = '../final-outputs/'

# ── Figure 1: Responsive rate by cancer type ─────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
cancer_rates = df.groupby('cancer_type')['predicted_responsive'].mean() * 100
colors = ['#1baf7a' if c == 'LIHC' else '#378add' if c == 'TNBC' else '#e24b4a' if c == 'PAAD' else '#ba7517' for c in cancer_rates.index]
bars = ax.bar(cancer_rates.index, cancer_rates.values, color=colors, width=0.5)
ax.set_ylabel('Predicted responsive rate (%)', fontsize=12)
ax.set_title('MENR predicted responsiveness by cancer type', fontsize=13)
ax.set_ylim(0, 35)
for bar, val in zip(bars, cancer_rates.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, f'{val:.1f}%', ha='center', fontsize=11)
plt.tight_layout()
plt.savefig(outdir + 'fig1_responsiveness_by_cancer.png', dpi=150)
plt.close()
print("Figure 1 saved.")

# ── Figure 2: Responsive rate by gene ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
gene_rates = df.groupby('gene')['predicted_responsive'].mean() * 100
colors = ['#1baf7a' if v > 20 else '#e24b4a' if v == 0 else '#378add' for v in gene_rates.values]
bars = ax.bar(gene_rates.index, gene_rates.values, color=colors, width=0.5)
ax.set_ylabel('Predicted responsive rate (%)', fontsize=12)
ax.set_title('MENR predicted responsiveness by gene', fontsize=13)
ax.set_ylim(0, 60)
for bar, val in zip(bars, gene_rates.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, f'{val:.1f}%', ha='center', fontsize=10)
plt.tight_layout()
plt.savefig(outdir + 'fig2_responsiveness_by_gene.png', dpi=150)
plt.close()
print("Figure 2 saved.")

# ── Figure 3: Methylation depth PAAD vs LIHC RASSF1 ─────────────────────────
paad = pd.read_csv('../data-acquisition/paad_rassf1_methylation_expression.tsv', sep='\t')
lihc = pd.read_csv('../data-acquisition/lihc_rassf1_methylation_expression.tsv', sep='\t')
mc = [c for c in paad.columns if 'Methylation' in c][0]
paad_heavy = paad[mc].dropna()
paad_heavy = paad_heavy[paad_heavy >= 0.3]
lihc_heavy = lihc[mc].dropna()
lihc_heavy = lihc_heavy[lihc_heavy >= 0.3]

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(paad_heavy, bins=20, alpha=0.6, color='#e24b4a', label=f'PAAD (n={len(paad_heavy)}, median={paad_heavy.median():.2f})')
ax.hist(lihc_heavy, bins=20, alpha=0.6, color='#1baf7a', label=f'LIHC (n={len(lihc_heavy)}, median={lihc_heavy.median():.2f})')
ax.set_xlabel('Starting RASSF1 methylation (beta value)', fontsize=12)
ax.set_ylabel('Number of patients', fontsize=12)
ax.set_title('RASSF1 methylation depth: PAAD vs LIHC\n(heavily silenced patients only, ≥0.3)', fontsize=12)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(outdir + 'fig3_paad_vs_lihc_rassf1_distribution.png', dpi=150)
plt.close()
print("Figure 3 saved.")

# ── Figure 4: Model B trajectory examples ────────────────────────────────────
def ode(t, m, kme=0.15, kde=0.40):
    dnmt = 0.12 * np.exp(-0.4*t) + 0.12*(1-0.85)
    tet1 = 0.4*(1-np.exp(-0.3*(t-12)))*np.exp(-0.08*(t-12)) if t>12 else 0
    return [kme*(1-m[0])*dnmt - kde*m[0]*tet1]

fig, ax = plt.subplots(figsize=(9, 5))
t_eval = np.linspace(0, 48, 200)
for m0, color, label in [(0.85, '#1baf7a', 'High methylation (0.85) — responsive'),
                          (0.65, '#378add', 'Medium methylation (0.65) — responsive'),
                          (0.40, '#e24b4a', 'Low methylation (0.40) — non-responsive'),
                          (0.30, '#ba7517', 'Borderline (0.30) — non-responsive')]:
    sol = solve_ivp(ode, (0, 48), [m0], t_eval=t_eval)
    final = sol.y[0][-1]
    pct = (m0 - final)/m0*100
    ax.plot(sol.t, sol.y[0], color=color, linewidth=2,
            label=f'{label} ({pct:.0f}% drop)')

ax.axhline(0.3, color='black', linestyle=':', linewidth=1, label='Example threshold')
ax.axvline(12, color='gray', linestyle='--', linewidth=1, label='Wave 2 onset (12hr)')
ax.set_xlabel('Time (hours)', fontsize=12)
ax.set_ylabel('Promoter methylation', fontsize=12)
ax.set_title('Model B: Two-wave epigenetic reprogramming trajectories', fontsize=12)
ax.legend(fontsize=9, loc='upper right')
ax.set_ylim(0, 1)
plt.tight_layout()
plt.savefig(outdir + 'fig4_model_b_trajectories.png', dpi=150)
plt.close()
print("Figure 4 saved.")

# ── Figure 5: Failure mode — genetic vs epigenetic ───────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
mech_rates = df.groupby('mechanism')['predicted_responsive'].mean() * 100
ax.bar(['Epigenetic silencing\n(MENR can act)', 'Genetic deletion\n(MENR cannot act)'],
       [mech_rates.get('epigenetic', 0), mech_rates.get('genetic', 0)],
       color=['#1baf7a', '#e24b4a'], width=0.4)
ax.set_ylabel('Predicted responsive rate (%)', fontsize=12)
ax.set_title('MENR failure mode: genetic vs epigenetic silencing', fontsize=12)
ax.set_ylim(0, 35)
ax.text(0, mech_rates.get('epigenetic', 0) + 0.5, f"{mech_rates.get('epigenetic', 0):.1f}%", ha='center', fontsize=12)
ax.text(1, 1, "0.0%", ha='center', fontsize=12)
plt.tight_layout()
plt.savefig(outdir + 'fig5_failure_mode_genetic_vs_epigenetic.png', dpi=150)
plt.close()
print("Figure 5 saved.")

print("\nAll figures saved to final-outputs/")
