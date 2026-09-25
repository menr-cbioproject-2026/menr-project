import pandas as pd

path = "../data-acquisition/GSE54503_RASSF1.csv"

try:
    df = pd.read_csv(path)
except Exception:
    df = pd.read_csv(path, sep='\t')

print("Columns found:", list(df.columns))

candidates = [c for c in df.columns if 'methyl' in c.lower() or 'beta' in c.lower()]
print("Likely methylation column(s):", candidates)

if candidates:
    col = candidates[0]
    values = df[col].dropna()
    print(f"\nUsing column: {col}")
    print(f"n = {len(values)}")
    print(f"Mean:   {values.mean():.4f}")
    print(f"Median: {values.median():.4f}")
else:
    print("\nNo obvious methylation column found — check the printed column list above and set 'col' manually.")
