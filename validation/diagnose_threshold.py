import numpy as np
from scipy.integrate import solve_ivp

def run_single(m0, kme, kde):
    def ode(t, m):
        dnmt = kme * np.exp(-0.1*t) + 0.05
        tet1 = kde * (1 - np.exp(-0.15*(t-12))) * np.exp(-0.05*t) if t > 12 else 0
        return -(kde*m[0]*tet1) + kme*(1-m[0])*dnmt
    sol = solve_ivp(ode, (0,96), [m0], t_eval=[96])
    final = sol.y[0][-1]
    pct_drop = (m0 - final)/m0 * 100
    return final, pct_drop

print('Starting m0   Final    % drop   Success (>=40% drop)')
for m0 in [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85]:
    final, drop = run_single(m0, 0.15, 0.4)
    print(f'{m0:.2f}         {final:.3f}    {drop:.1f}%    {drop >= 40}')
