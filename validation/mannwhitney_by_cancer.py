import pandas as pd
from scipy import stats

df = pd.read_csv('../digital-twin/digital_twin_results_v2.csv')
epi = df[df['mechanism']=='epigenetic']

print(f"{'Cancer':<8} {'Resp median':<14} {'Non-resp median':<18} {'p-value':<12} {'Significant?'}")
print("-"*65)

for cancer in ['LIHC','TNBC','LUAD','PAAD']:
    sub = epi[epi['cancer_type']==cancer]
    resp = sub[sub['predicted_responsive']==True]['starting_methylation'].dropna()
    nonresp = sub[sub['predicted_responsive']==False]['starting_methylation'].dropna()
    if len(resp) > 0 and len(nonresp) > 0:
        stat, p = stats.mannwhitneyu(resp, nonresp, alternative='greater')
        sig = 'YES ***' if p < 0.001 else ('YES *' if p < 0.05 else 'NO')
        print(f"{cancer:<8} {resp.median():<14.3f} {nonresp.median():<18.3f} {p:<12.2e} {sig}")
