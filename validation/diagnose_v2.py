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

# Test the hardest combo (Try 4) on the WORST case patient (0.95 start)
t_eval = np.linspace(0, 96, 20)
sol = solve_ivp(methylation_ode, (0, 96), [0.95], args=(0.35, 0.08), t_eval=t_eval, method='RK45')

print("Worst-case patient (start=0.95), hardest settings (kme=0.35, kde=0.08):")
print("Time(hr)  Methylation")
for t, m in zip(sol.t, sol.y[0]):
    print(f"{t:6.1f}    {m:.4f}")

# What is the LONG-TERM resting point (steady state) once TET1 clears out?
# After TET1 fades away (~hour 60+), what does methylation settle at?
print(f"\nFinal value at hour 96: {sol.y[0][-1]:.4f}")
print(f"Reactivation threshold: 0.30")
print(f"Did it cross? {'YES' if sol.y[0][-1] < 0.3 else 'NO'}")
