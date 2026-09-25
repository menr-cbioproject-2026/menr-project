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

# Test multiple starting points at the DOCUMENT-CORRECT 48-hour window
test_starts = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95]
t_eval = np.linspace(0, 48, 10)

print("Using the project document's actual stated Wave 2 window (0-48hrs), kme=0.25, kde=0.18:\n")
print(f"{'Start':<8} {'Final @ 48hr':<14} {'Crosses 0.3?'}")
for start in test_starts:
    sol = solve_ivp(methylation_ode, (0, 48), [start], args=(0.25, 0.18), t_eval=t_eval, method='RK45')
    final = sol.y[0][-1]
    print(f"{start:<8} {final:<14.4f} {'YES' if final < 0.3 else 'NO'}")
