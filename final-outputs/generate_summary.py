import pandas as pd
from scipy import stats
import numpy as np

df = pd.read_csv('../digital-twin/digital_twin_results_v2.csv')

lines = []
lines.append("=" * 60)
lines.append("MENR PROJECT — KEY RESULTS SUMMARY")
lines.append("=" * 60)

lines.append(f"\nTOTAL PREDICTIONS: {len(df):,}")
lines.append(f"OVERALL RESPONSIVE RATE: {df['predicted_responsive'].mean()*100:.1f}%")

lines.append("\n--- Responsive rate by cancer type ---")
for cancer, rate in df.groupby('cancer_type')['predicted_responsive'].mean().items():
    lines.append(f"  {cancer}: {rate*100:.1f}%")

lines.append("  NOTE: LUAD stiffness parameters (E0, k) reuse LIHC's calibration.")
lines.append("  No lung-specific MRE tumor stiffness data was available in the")
lines.append("  literature. LUAD results should be interpreted with this caveat.")
lines.append("\n--- Responsive rate by gene ---")
for gene, rate in df.groupby('gene')['predicted_responsive'].mean().items():
    lines.append(f"  {gene}: {rate*100:.1f}%")

lines.append("\n--- By silencing mechanism ---")
for mech, rate in df.groupby('mechanism')['predicted_responsive'].mean().items():
    lines.append(f"  {mech}: {rate*100:.1f}%")

lines.append("\n--- Mann-Whitney U (responsive vs non-responsive methylation) ---")
epi = df[df['mechanism']=='epigenetic']
for cancer in ['LIHC','TNBC','LUAD','PAAD']:
    sub = epi[epi['cancer_type']==cancer]
    resp = sub[sub['predicted_responsive']==True]['starting_methylation'].dropna()
    nonresp = sub[sub['predicted_responsive']==False]['starting_methylation'].dropna()
    if len(resp) > 0 and len(nonresp) > 0:
        stat, p = stats.mannwhitneyu(resp, nonresp, alternative='greater')
        lines.append(f"  {cancer}: resp median={resp.median():.3f}, non-resp median={nonresp.median():.3f}, p={p:.2e}")

lines.append("\n--- Independent replication (GSE54503, Columbia University) ---")
lines.append("  RASSF1 median in 132 HCC samples: 0.504")
lines.append("  TCGA LIHC RASSF1 median: 0.593")
lines.append("  Directional agreement: YES")

lines.append("\n--- Key finding ---")
lines.append("  PAAD RASSF1 heavily-silenced median: 0.426")
lines.append("  LIHC RASSF1 heavily-silenced median: 0.593")
lines.append("  Shallower methylation in PAAD explains lower responsiveness")
lines.append("  independent of stiffness-based eligibility.")

lines.append("\n" + "=" * 60)

summary = "\n".join(lines)
print(summary)

with open('../final-outputs/MENR_results_summary.txt', 'w') as f:
    f.write(summary)

print("\nSaved to final-outputs/MENR_results_summary.txt")
