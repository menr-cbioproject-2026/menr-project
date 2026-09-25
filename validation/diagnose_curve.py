import numpy as np
from scipy.integrate import solve_ivp

def DNMT3A_active(t, decay=0.3):
    return np.exp(-decay * t)

def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

def methylation_ode(t, m, kme=0.15, kde=0.4):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]

# Test with a HIGH starting methylation (worst case)
t_eval = np.linspace(0, 96, 20)
sol = solve_ivp(methylation_ode, (0, 96), [0.95], t_eval=t_eval, method='RK45')

print("Time(hr)  Methylation")
for t, m in zip(sol.t, sol.y[0]):
    print(f"{t:6.1f}    {m:.4f}")
