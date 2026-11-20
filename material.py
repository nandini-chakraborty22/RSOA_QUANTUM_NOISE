"""
Bulk InGaAsP material-gain model (parabolic bands), following
M. J. Connelly, IEEE J. Quantum Electron. 37(3), 439-447 (2001), Eqs. (2)-(17).

All quantities SI unless stated. Energies inside are in joules.
"""
import numpy as np

# Physical constants
q = 1.602176634e-19
h = 6.62607015e-34
hbar = h / (2 * np.pi)
c = 299792458.0
kB = 1.380649e-23

# ---------------------------------------------------------------------------
# Material / geometry parameters (bulk InP-InGaAsP buried-ridge SOA, after Connelly 2001 Table I).
# Users should verify against the original table before publication.
P = dict(
    y=0.892,            # As molar fraction
    a_eg=1.35, b_eg=-0.775, c_eg=0.149,   # Eg0 = a + b*y + c*y^2  [eV]
    Kg=0.9e-10,         # bandgap shrinkage coefficient [eV m]
    me=4.10e-32,        # CB electron effective mass [kg]
    mhh=4.19e-31,       # VB heavy-hole effective mass [kg]
    mlh=5.06e-32,       # VB light-hole effective mass [kg]
    n1=3.22,            # active-region refractive index
    T=300.0,            # temperature [K]
    Arad=1e7, Brad=5.6e-16,               # radiative recombination
    Anrad=3.5e8, Bnrad=0.0, Caug=3e-41,   # non-radiative recombination
    K0=6200.0, K1=7500e-24,               # material loss  alpha = K0 + Gamma*K1*n
    Gamma=0.45,         # optical confinement factor
    W=0.4e-6, d=0.4e-6, # active width and thickness [m]
)


def Eg(n, p=P):
    """Bandgap [J] including carrier-induced shrinkage."""
    y = p['y']
    eg0 = p['a_eg'] + p['b_eg'] * y + p['c_eg'] * y ** 2
    return q * (eg0 - p['Kg'] * np.cbrt(n))


def _nilsson(x):
    return np.log(x) + x * (64 + 0.05524 * x * (64 + np.sqrt(x))) ** (-0.25)


def quasi_fermi(n, p=P):
    kT = kB * p['T']
    nc = 2 * (p['me'] * kT / (2 * np.pi * hbar ** 2)) ** 1.5
    mdh = (p['mhh'] ** 1.5 + p['mlh'] ** 1.5) ** (2 / 3)
    nv = 2 * (mdh * kT / (2 * np.pi * hbar ** 2)) ** 1.5
    Efc = _nilsson(n / nc) * kT          # above CB edge
    Efv = -_nilsson(n / nv) * kT         # relative to VB top (negative = into band)
    return Efc, Efv


def tau_rad(n, p=P):
    """Radiative carrier lifetime (Connelly Eq. 52)."""
    return 1.0 / (p['Arad'] + p['Brad'] * n)


def gain_components(nu, n, p=P):
    """
    Return (g_m, g_e) [1/m]: net material gain and the emission component used
    for spontaneous-emission coupling.  nu [Hz] and n [m^-3] broadcast.
    Eq. (1) of the revised manuscript:
        g_m = c^2/(4*sqrt(2)*pi^(3/2)*n1^2*tau*nu^2) * (2 me mhh/(hbar(me+mhh)))^(3/2)
              * sqrt(nu - Eg/h) * (fc - fv)
    """
    nu = np.asarray(nu, float)
    n = np.asarray(n, float)
    me, mhh, kT = p['me'], p['mhh'], kB * p['T']
    eg = Eg(n, p)
    Efc, Efv = quasi_fermi(n, p)
    dE = h * nu - eg
    above = dE > 0
    dEp = np.where(above, dE, 0.0)
    Ea = dEp * mhh / (me + mhh)
    Eb = -dEp * me / (me + mhh)
    fc = 1.0 / (1.0 + np.exp((Ea - Efc) / kT))
    fv = 1.0 / (1.0 + np.exp((Eb - Efv) / kT))
    pref = c ** 2 / (4 * np.sqrt(2) * np.pi ** 1.5 * p['n1'] ** 2 * tau_rad(n, p) * nu ** 2)
    pref = pref * (2 * me * mhh / (hbar * (me + mhh))) ** 1.5 * np.sqrt(dEp / h)
    gm = np.where(above, pref * (fc - fv), 0.0)
    ge = np.where(above, pref * fc * (1 - fv), 0.0)
    return gm, ge


def recombination(n, p=P):
    return ((p['Arad'] + p['Anrad']) * n + (p['Brad'] + p['Bnrad']) * n ** 2
            + p['Caug'] * n ** 3)


def dR_dn(n, p=P):
    return ((p['Arad'] + p['Anrad']) + 2 * (p['Brad'] + p['Bnrad']) * n
            + 3 * p['Caug'] * n ** 2)


def alpha_loss(n, p=P):
    return p['K0'] + p['Gamma'] * p['K1'] * n
