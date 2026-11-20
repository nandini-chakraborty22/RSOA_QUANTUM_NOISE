"""Noise-figure tools for the steady-state RSOA model (rsoa_steady.py).

The chip-internal noise figure is computed in two ways:
  route A  distributed: the double pass is unfolded into 2K sections and the
           ASE born in each section is propagated to the output with the net
           gain that follows it (exact per-section solution, so sections with
           net loss are handled without the n_sp*g/(g - alpha) singularity);
  route B  from the simulated ASE: NF = 2*S/G + 1/G, with S the output ASE
           spectral density in the signal polarisation in units of h*nu.
Chip-internal NF excludes fibre coupling loss, which adds to the
fibre-to-fibre NF in dB but is not amplifier physics.

A narrow probe bin is added to the ASE grid exactly at the signal frequency
so that route B needs no interpolation.
"""
import numpy as np
import material as m
from rsoa_steady import RSOA

PROBE_BW = 1e9      # Hz; narrow enough not to perturb the carrier density
QUANTUM_LIMIT_DB = 10 * np.log10(2)


def db(x):
    return 10 * np.log10(x)


def dbm2w(p_dbm):
    return 1e-3 * 10 ** (np.asarray(p_dbm, float) / 10)


def make_device(L=600e-6, R1=1e-5, R2=0.9, I=0.08, lam=1550e-9, n_ase=48):
    """Reference RSOA with a probe bin at the signal wavelength."""
    d = RSOA(L=L, R1=R1, R2=R2, I=I, lam_s=lam, n_ase=n_ase)
    d.nu_a = np.append(d.nu_a, d.nu_s)
    d.dnu = np.append(d.dnu, PROBE_BW)
    d.j_sig = len(d.nu_a) - 1
    return d


def solve(d, p_dbm, n0=None):
    s = d.solve(dbm2w(p_dbm), n_init=n0, tol=1e-5)
    s['Pin'] = float(dbm2w(p_dbm))
    return s


def nf_distributed(dz, ggm, gge, alpha, R1, R2):
    """Route A. ggm, gge, alpha: per-section Gamma*g_m, Gamma*g_e and alpha at the
    signal (1/m), section 0 at the front facet. Returns (NF, G), chip-internal, linear."""
    K = len(ggm)
    path = list(range(K)) + ['rear'] + list(range(K - 1, -1, -1))
    S, G = 0.0, 1.0
    for k in path:
        if k == 'rear':
            S *= R2; G *= R2; continue
        g = ggm[k] - alpha[k]
        e = np.exp(g * dz)
        S = S * e + (gge[k] * (e - 1) / g if abs(g * dz) > 1e-12 else gge[k] * dz)
        G *= e
    S *= (1 - R1); G *= (1 - R1)          # front-facet transmission on the way out
    return 2 * S / G + 1 / G, G


def net_transparency_density(lam=1550e-9, p=m.P):
    """Carrier density at which modal gain equals internal loss at lam."""
    nu = m.c / lam
    lo, hi = 1.0e24, 3.0e24
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        gm, _ = m.gain_components(nu, mid, p)
        if p['Gamma'] * gm > m.alpha_loss(mid, p): hi = mid
        else: lo = mid
    return mid


def analyse(d, s):
    """Both NF routes plus diagnostics for one solved operating point."""
    p, hv = d.p, m.h * d.nu_s
    Qin = d.eta * (1 - d.R1) * s['Pin'] / hv               # photons/s inside the chip
    G_chip = (1 - d.R1) * s['Sb'][0] / Qin
    j = d.j_sig
    Ab0 = s['ase_bins'][j] / (2 * d.eta * (1 - d.R1) * hv)  # photons/s, one polarisation
    S_out = (1 - d.R1) * Ab0 / d.dnu[j]                    # photons/s/Hz leaving the chip
    nfB = 2 * S_out / G_chip + 1 / G_chip
    gm, ge = m.gain_components(d.nu_s, d.n, p)
    Gam = p['Gamma']
    alpha = m.alpha_loss(d.n, p)
    nfA, GA = nf_distributed(d.dz, Gam * gm, Gam * ge, alpha, d.R1, d.R2)
    return dict(G_chip_dB=float(db(G_chip)), G_routeA_dB=float(db(GA)),
                G_fibre_dB=float(db(s['G_ff'])),
                NF_A_dB=float(db(nfA)), NF_B_dB=float(db(nfB)),
                NF_fibre_dB=float(db(nfA) - db(d.eta)),
                nsp_eff_gain_weighted=float((nfA - 1 / GA) * GA / (2 * (GA - 1))),
                sections_net_loss=int(np.sum(Gam * gm < alpha)),
                n_front6_mean=float(d.n[:6].mean()), n_mean=float(d.n.mean()))
