with open('digital_twin.py', 'r') as f:
    content = f.read()

old = """def DNMT3A_active(t, natural_level=0.05, inhibitor_strength=0.85, inhibitor_decay=0.12):
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)
def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))
def methylation_ode(t, m, kme=0.20, kde=0.30):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]
RECOVERY_FRACTION = 0.40
def model_b_outcome(start_m):
    t_span = (0, 96)
    t_eval = np.linspace(0, 96, 100)
    solution = solve_ivp(methylation_ode, t_span, [start_m], t_eval=t_eval, method='RK45')"""

new = """def DNMT3A_active(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)
def TET1_active(t, delay=12, rise=0.3, clearance=0.08):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))
def methylation_ode(t, m, kme=0.15, kde=0.40):
    dnmt3a = DNMT3A_active(t)
    tet1   = TET1_active(t)
    return [kme * (1 - m) * dnmt3a - kde * m * tet1]
RECOVERY_FRACTION = 0.40
def model_b_outcome(start_m):
    t_span = (0, 48)
    t_eval = np.linspace(0, 48, 100)
    solution = solve_ivp(methylation_ode, t_span, [start_m], t_eval=t_eval, method='RK45')"""

if old in content:
    content = content.replace(old, new)
    with open('digital_twin.py', 'w') as f:
        f.write(content)
    print("Patched successfully.")
else:
    print("ERROR: old string not found - check spacing")
