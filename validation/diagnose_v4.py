import numpy as np
from scipy.integrate import solve_ivp

def DNMT3A_active(t, decay=0.3):
    return np.exp(-decay * t)

def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

def methylation_ode(t, m, kme, kde):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]

test_starts = [0.3, 0.5, 0.7, 0.95]
t_eval = np.linspace(0, 96, 10)

# Try MUCH more balanced rates - kme and kde close to equal
combos = [
    (0.4, 0.05, "kme=0.4, kde=0.05 (DNMT3A much stronger)"),
    (0.3, 0.1,  "kme=0.3, kde=0.1"),
    (0.25, 0.15,"kme=0.25, kde=0.15 (closer balance)"),
]

for kme, kde, label in combos:
    print(f"\n--- {label} ---")
    print(f"{'Start':<8} {'Final @ 96hr':<14} {'Crosses 0.3?'}")
    for start in test_starts:
        sol = solve_ivp(methylation_ode, (0, 96), [start], args=(kme, kde), t_eval=t_eval, method='RK45')
        final = sol.y[0][-1]
        print(f"{start:<8} {final:<14.4f} {'YES' if final < 0.3 else 'NO'}")
