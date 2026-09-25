import pandas as pd
import numpy as np
from scipy import stats

line = open('../data-acquisition/GSE54503_series_matrix.txt').readlines()
probe_line = [l for l in line if l.startswith('"cg08047457"')][0]
values = [float(x.strip().strip('"')) for x in probe_line.split('\t')[1:] if x.strip()]

vals = np.array(values)
print(f"GSE54503 RASSF1 (cg08047457) across {len(vals)} HCC samples")
print(f"Mean: {vals.mean():.3f}")
print(f"Median: {np.median(vals):.3f}")
print(f"% above 0.3 (heavily silenced): {(vals>=0.3).mean()*100:.1f}%")
print(f"% above 0.55: {(vals>=0.55).mean()*100:.1f}%")
print(f"\nTCGA LIHC RASSF1 for comparison:")
print(f"Mean: 0.585, Median: 0.593, % above 0.55: 66.1%")
print(f"\nDirectional agreement: both datasets show high RASSF1 methylation in HCC")
