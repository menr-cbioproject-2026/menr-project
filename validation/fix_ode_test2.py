import numpy as np
from scipy.integrate import solve_ivp

def run_single(m0, kme=0.15, kde=0.4):
    def ode(t, m):
        dnmt_active = 0.3 * np.exp(-0.5*t) + 0.08
        tet1 = 0.6 * (t/12) * np.exp(-0.04*(t-12)) if t > 12 else 0
        dmdt = kme*(1-m[0])*dnmt_active - kde*m[0]*tet1
        return dmdt
    sol = solve_ivp(ode, (0,96), [m0], t_eval=[96], max_step=0.5)
    final = max(0, min(1, sol.y[0][-1]))
    pct_drop = (m0 - final)/m0 * 100
    return final, pct_drop

print('Starting m0   Final    % drop   Success (>=40% drop)')
for m0 in [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85, 0.95]:
    final, drop = run_single(m0)
    print(f'{m0:.2f}         {final:.3f}    {drop:.1f}%    {drop >= 40}')
