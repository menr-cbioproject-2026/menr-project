with open('../model-b-kinetics/model_b.py', 'r') as f:
    content = f.read()

old_dnmt = """def DNMT3A_active(t, natural_level=0.05, inhibitor_strength=0.85, inhibitor_decay=0.12):
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)"""

new_dnmt = """def DNMT3A_active(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
    # natural_level raised to 0.12 (from 0.05) so background DNMT3A provides
    # meaningful competition against TET1 - tuned fit parameter, stated uncertainty
    # inhibitor_decay raised to 0.4 so Wave 1 suppression fades realistically
    inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
    return natural_level * (1 - inhibition)"""

old_tet1 = """def TET1_active(t, delay=12, rise=0.3, clearance=0.05):
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))"""

new_tet1 = """def TET1_active(t, delay=12, rise=0.3, clearance=0.08):
    # clearance raised to 0.08 (from 0.05) so TET1 fades within 48hr window
    # this creates the biologically realistic gradient: high starters succeed
    # before TET1 clears, low starters do not - tuned fit parameter
    if t < delay:
        return 0.0
    return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))"""

old_ode = """def methylation_ode(t, m, kme=0.20, kde=0.30):"""
new_ode = """def methylation_ode(t, m, kme=0.15, kde=0.40):"""

old_window = """solve_ivp(methylation_ode, (0, 96)"""
new_window = """solve_ivp(methylation_ode, (0, 48)"""

content = content.replace(old_dnmt, new_dnmt)
content = content.replace(old_tet1, new_tet1)
content = content.replace(old_ode, new_ode)
content = content.replace(old_window, new_window)

with open('../model-b-kinetics/model_b.py', 'w') as f:
    f.write(content)

print("Done. Verify changes:")
with open('../model-b-kinetics/model_b.py', 'r') as f:
    lines = f.readlines()
for i, line in enumerate(lines[23:55], start=23):
    print(f"{i}: {line}", end='')
