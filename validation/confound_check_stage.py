import pandas as pd

STAGE_COL = 'Neoplasm Disease Stage American Joint Committee on Cancer Code'

def load_clinical(path):
    df = pd.read_csv(path, sep='\t')
    df = df[['Sample ID', STAGE_COL]].rename(columns={STAGE_COL: 'stage_raw'})
    return df

def load_methylation(path):
    df = pd.read_csv(path, sep='\t')
    meth_col = [c for c in df.columns if 'Methylation' in c][0]
    df = df[['Sample Id', meth_col]].rename(columns={'Sample Id': 'Sample ID', meth_col: 'methylation'})
    df = df.dropna(subset=['methylation'])
    return df

def simplify_stage(s):
    if pd.isna(s):
        return None
    s = str(s).upper()
    if 'IV' in s:
        return 'IV'
    if 'III' in s:
        return 'III'
    if 'II' in s:
        return 'II'
    if 'I' in s:
        return 'I'
    return None

results = {}
for cancer, clin_path, meth_path in [
    ('LIHC', 'data-acquisition/lihc_tcga_pan_can_atlas_2018_clinical_data.tsv',
             'data-acquisition/lihc_rassf1_methylation_expression.tsv'),
    ('PAAD', 'data-acquisition/paad_tcga_pan_can_atlas_2018_clinical_data.tsv',
             'data-acquisition/paad_rassf1_methylation_expression.tsv'),
]:
    clin = load_clinical(clin_path)
    meth = load_methylation(meth_path)
    merged = pd.merge(meth, clin, on='Sample ID')
    merged['stage_simple'] = merged['stage_raw'].apply(simplify_stage)
    merged = merged.dropna(subset=['stage_simple'])

    print(f"\n=== {cancer} ===")
    print(f"Methylation rows: {len(meth)} | Matched+staged: {len(merged)}")
    print(f"Stage distribution:\n{merged['stage_simple'].value_counts().sort_index()}")
    print(f"Median methylation by stage:")
    print(merged.groupby('stage_simple')['methylation'].agg(['median', 'count']))

    results[cancer] = merged

print("\n=== COMPARISON: LIHC vs PAAD within matching stages ===")
for stage in ['I', 'II', 'III', 'IV']:
    lihc_sub = results['LIHC'][results['LIHC']['stage_simple'] == stage]['methylation']
    paad_sub = results['PAAD'][results['PAAD']['stage_simple'] == stage]['methylation']
    if len(lihc_sub) > 0 and len(paad_sub) > 0:
        print(f"Stage {stage}: LIHC median={lihc_sub.median():.3f} (n={len(lihc_sub)}), "
              f"PAAD median={paad_sub.median():.3f} (n={len(paad_sub)})")
    else:
        print(f"Stage {stage}: insufficient data (LIHC n={len(lihc_sub)}, PAAD n={len(paad_sub)})")
