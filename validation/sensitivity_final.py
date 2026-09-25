import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

df = pd.read_csv('../digital-twin/digital_twin_results_v2.csv')
borderline = df[(df['mechanism']=='epigenetic') & 
                (df['starting_methylation'] >= 0.3) & 
                (df['starting_methylation'] <= 0.55)]['starting_methylation'].dropna().values

def run_model(kme, kde, starts):
    results = []
    for m0 in starts:
        def ode(t, m):
            dnmt = 0.15 * np.exp(-0.1*t) + 0.05
            tet1 = 0.4 * (1 - np.exp(-0.15*(t-12))) * np.exp(-0.05*t) if t > 12 else 0
            return kme*(1-m[0])*dnmt - kde*m[0]*tet1
        sol = solve_ivp(ode, (0,96), [m0], t_eval=[96])
        final = sol.y[0][-1]
        results.append(final < m0 * 0.6)
    return np.mean(results) * 100

baseline = run_model(0.15, 0.4, borderline)
print(f"Borderline patients (0.3-0.55 methylation): n={len(borderline)}")
print(f"\n{'Scenario':<30} {'% responsive':<15} {'Shift'}")
print("-"*55)
print(f"{'Baseline (kme=0.15, kde=0.4)':<30} {baseline:<15.1f} --")

scenarios = [
    ("kme x2 (DNMT3A stronger)", 0.30, 0.4),
    ("kme x3 (DNMT3A much stronger)", 0.45, 0.4),
    ("kde x0.5 (TET1 weaker)", 0.15, 0.2),
    ("kde x0.25 (TET1 much weaker)", 0.15, 0.1),
    ("both adverse", 0.30, 0.2),
    ("both favorable", 0.10, 0.6),
]

for label, kme, kde in scenarios:
    result = run_model(kme, kde, borderline)
    shift = result - baseline
    print(f"{label:<30} {result:<15.1f} {shift:+.1f} pts")
