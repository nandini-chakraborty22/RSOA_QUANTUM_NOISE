"""Noise figure of the reference RSOA (600 um, 80 mA, R1 = 1e-5, R2 = 0.9, 1550 nm)
compared with the 3 dB quantum limit.
Writes results/noise_figure.json and figures/FigNF_noise_figure.(png|pdf)."""
import json, os, time
import numpy as np
import material as m
from rsoa_steady import RSOA
import noise as nz
import style
from style import plt, C, panel, save

os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
style.setup()
t0 = time.time()
res = {}

# 1. NF versus input power
d = nz.make_device()
Pins = np.arange(-40.0, 3.01, 1.0)
sweep, n0 = [], None
for P in Pins:
    s = nz.solve(d, P, n0); n0 = s['n'].copy()
    r = nz.analyse(d, s); r['Pin_dBm'] = float(P); sweep.append(r)
res['sweep'] = sweep
Gf = np.array([r['G_fibre_dB'] for r in sweep])
P3 = float(np.interp(Gf[0] - 3, Gf[::-1], Pins[::-1]))
key = {}
for name, P in [('linear', -40.0), ('compression_3dB', P3), ('seed', -13.0), ('deep_saturation', 0.0)]:
    r = nz.analyse(d, nz.solve(d, P)); r['Pin_dBm'] = float(P); key[name] = r
res['key_points'] = key
res['max_route_difference_dB'] = float(max(abs(r['NF_A_dB'] - r['NF_B_dB']) for r in sweep))

# effect of the probe bin on the gain of the unmodified model
d0 = RSOA(L=600e-6, R1=1e-5, R2=0.9, I=0.08)
res['probe_bin_gain_change_dB'] = float(key['linear']['G_fibre_dB']
                                        - nz.db(d0.solve(nz.dbm2w(-40.0), tol=1e-5)['G_ff']))

# 2. Decomposition of the gap at -40 dBm
p = m.P
Jvol = d.I / (m.q * d.V)
lo, hi = 1e22, 1e26
for _ in range(200):
    mid = np.sqrt(lo * hi)
    if m.recombination(mid, p) > Jvol: hi = mid
    else: lo = mid
n_bias = mid                                  # carrier density with no light at all
gm0, ge0 = m.gain_components(d.nu_s, n_bias, p)
K = d.K
ggm = np.full(K, p['Gamma'] * gm0); gge = np.full(K, p['Gamma'] * ge0)
nfI, _ = nz.nf_distributed(d.dz, ggm, gge, np.zeros(K), d.R1, d.R2)
nfL, _ = nz.nf_distributed(d.dz, ggm, gge, np.full(K, m.alpha_loss(n_bias, p)), d.R1, d.R2)
res['decomposition'] = dict(n_bias=float(n_bias), nsp_bias=float(ge0 / gm0),
                            ideal_no_loss_NF_dB=float(nz.db(nfI)),
                            with_loss_NF_dB=float(nz.db(nfL)),
                            full_model_NF_dB=key['linear']['NF_A_dB'])
res['n_net_transparency_1550'] = float(nz.net_transparency_density())
res['n_material_transparency_1550'] = 1.15e24

# carrier profiles
prof = {}
for P in (-40, -13, 0):
    nz.solve(d, P); prof[str(P)] = d.n.tolist()
res['profiles_n'] = prof
json.dump(res, open('results/noise_figure.json', 'w'), indent=1)

# Figure
fig, ax = plt.subplots(1, 3, figsize=(style.COL2, 2.3))
P = [r['Pin_dBm'] for r in sweep]
ax[0].plot(P, [r['NF_A_dB'] for r in sweep], color=C[0], label='Route A (distributed)')
ax[0].plot(P[::3], [r['NF_B_dB'] for r in sweep][::3], 'o', color=C[1], mfc='none', label='Route B (simulated ASE)')
ax[0].axhline(nz.QUANTUM_LIMIT_DB, color='k', ls='--', lw=0.8, label='Quantum limit (3 dB)')
ax[0].axvline(-13, color=C[5], ls=':', lw=0.8)
ax[0].set_xlabel('Input power (dBm)'); ax[0].set_ylabel('Chip-internal NF (dB)')
ax[0].set_ylim(0, 16); ax[0].legend(loc='center left', bbox_to_anchor=(0.0, 0.55)); panel(ax[0], 'a')

sec = np.arange(1, K + 1)
for (Pk, col) in [('-40', C[0]), ('-13', C[3]), ('0', C[1])]:
    ax[1].plot(sec, np.array(prof[Pk]) / 1e24, color=col, label=f'{Pk} dBm')
ax[1].axhline(res['n_net_transparency_1550'] / 1e24, color='k', ls=':', lw=0.8,
              label='net transparency, 1550 nm')
ax[1].set_xlabel('Section (1 = front facet)'); ax[1].set_ylabel(r'$n$ (10$^{24}$ m$^{-3}$)')
ax[1].legend(fontsize=6); panel(ax[1], 'b')

dec = res['decomposition']
vals = [dec['ideal_no_loss_NF_dB'], dec['with_loss_NF_dB'], dec['full_model_NF_dB']]
ax[2].bar(['Ideal\ninversion', '+ internal\nloss', 'Full model\n(+ ASE)'], vals,
          color=[C[2], C[3], C[0]], width=0.6)
ax[2].axhline(nz.QUANTUM_LIMIT_DB, color='k', ls='--', lw=0.8)
for i, v in enumerate(vals):
    ax[2].text(i, v + 0.3, f'{v:.1f}', ha='center', fontsize=7)
ax[2].set_ylabel('NF at -40 dBm (dB)'); ax[2].set_ylim(0, 16); ax[2].grid(axis='x')
panel(ax[2], 'c')
fig.tight_layout(); save(fig, 'FigNF_noise_figure', 'figures')

print(f'done in {time.time() - t0:.1f} s')
for k, v in key.items():
    print(f"{k:16s} Pin {v['Pin_dBm']:6.2f} dBm  NF {v['NF_A_dB']:.2f} / {v['NF_B_dB']:.2f} dB  "
          f"net-loss sections {v['sections_net_loss']}")
print('decomposition', {k: round(v, 3) for k, v in dec.items()})
print('max route difference (dB)', res['max_route_difference_dB'])
print('probe-bin gain change (dB)', res['probe_bin_gain_change_dB'])
