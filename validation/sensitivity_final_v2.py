import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

df = pd.read_csv('../digital-twin/digital_twin_results_v2.csv')
borderline = df[
    (df['mechanism']=='epigenetic') & 
    (df['starting_methylation'] >= 0.3) & 
    (df['starting_methylation'] <= 0.65)
]['starting_methylation'].dropna().values

def run_model(kme, kde, starts):
    results = []
    for m0 in starts:
        def ode(t, m, kme=kme, kde=kde):
            dnmt = 0.12 * np.exp(-0.4*t) + 0.12*(1-0.85)
            tet1 = 0.4*(1-np.exp(-0.3*(t-12)))*np.exp(-0.08*(t-12)) if t>12 else 0
            return [kme*(1-m[0])*dnmt - kde*m[0]*tet1]
        sol = solve_ivp(ode, (0,48), [m0], t_eval=[48], max_step=0.5)
        final = sol.y[0][-1]
        pct_drop = (m0-final)/m0*100
        results.append(pct_drop >= 40)
    return np.mean(results)*100

baseline = run_model(0.15, 0.40, borderline)
print(f"Borderline patients (0.3-0.65): n={len(borderline)}")
print(f"\n{'Scenario':<35} {'% responsive':<15} {'Shift'}")
print("-"*60)
print(f"{'Baseline (kme=0.15, kde=0.40)':<35} {baseline:<15.1f} --")

scenarios = [
    ("kme x2 — DNMT3A stronger",       0.30, 0.40),
    ("kme x3 — DNMT3A much stronger",   0.45, 0.40),
    ("kde x0.5 — TET1 weaker",          0.15, 0.20),
    ("kde x0.25 — TET1 much weaker",    0.15, 0.10),
    ("both adverse",                    0.30, 0.20),
    ("both favorable",                  0.10, 0.60),
]

for label, kme, kde in scenarios:
    result = run_model(kme, kde, borderline)
    shift = result - baseline
    print(f"{label:<35} {result:<15.1f} {shift:+.1f} pts")
