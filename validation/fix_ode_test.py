import numpy as np
from scipy.integrate import solve_ivp

def run_single(m0, kme=0.15, kde=0.4):
    def ode(t, m):
        dnmt_active = np.exp(-0.3*t) + 0.05
        tet1 = (1 - np.exp(-0.2*(t-12))) * np.exp(-0.03*t) if t > 12 else 0
        dmdt = kme*(1-m[0])*dnmt_active - kde*m[0]*tet1
        return dmdt
    sol = solve_ivp(ode, (0,96), [m0], t_eval=[96])
    final = sol.y[0][-1]
    pct_drop = (m0 - final)/m0 * 100
    return final, pct_drop

print('Starting m0   Final    % drop   Success (>=40% drop)')
for m0 in [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85]:
    final, drop = run_single(m0)
    print(f'{m0:.2f}         {final:.3f}    {drop:.1f}%    {drop >= 40}')
