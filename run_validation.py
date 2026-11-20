"""Validation checks for the noise-figure result. Writes results/validation.json and
figures/FigV_validation.(png|pdf).

(1) Rear-facet reflectivity sweep: if ASE returning from the rear facet clamps the
    front sections, lowering R2 should raise the front-section carrier density
    well above net transparency and remove the gap.
(2) ASE grid refinement: the result should not depend on the number of ASE bins.
(3) Signal-wavelength sweep, 1530-1570 nm: is the gap specific to 1550 nm?
All at -40 dBm input (linear regime), 600 um, 80 mA, R1 = 1e-5.
"""
import json, os, time
import numpy as np
import material as m
import noise as nz
import style
from style import plt, C, panel, save

os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
style.setup()
t0 = time.time()
P_IN = -40.0
res = {}

# (1) rear-facet reflectivity
rows, prof = [], {}
for R2 in [0.9, 0.5, 0.3, 0.1, 0.03, 0.01]:
    d = nz.make_device(R2=R2)
    s = nz.solve(d, P_IN)
    r = nz.analyse(d, s); r['R2'] = R2
    rows.append(r); prof[str(R2)] = d.n.tolist()
res['rear_facet'] = rows
res['rear_facet_profiles'] = prof
print('(1) done', time.time() - t0)

# (2) ASE grid refinement
rows = []
for nb in [24, 48, 96, 192]:
    d = nz.make_device(n_ase=nb)
    r = nz.analyse(d, nz.solve(d, P_IN)); r['n_ase_bins'] = nb
    rows.append(r)
res['ase_grid'] = rows
print('(2) done', time.time() - t0)

# (3) signal wavelength
rows = []
for lam_nm in [1530, 1540, 1550, 1560, 1570]:
    d = nz.make_device(lam=lam_nm * 1e-9)
    r = nz.analyse(d, nz.solve(d, P_IN)); r['lam_nm'] = lam_nm
    rows.append(r)
res['wavelength'] = rows
print('(3) done', time.time() - t0)

json.dump(res, open('results/validation.json', 'w'), indent=1)

# Figure
fig, axs = plt.subplots(2, 2, figsize=(style.COL2, 4.4))
ax = axs.ravel()
rf = res['rear_facet']
R2s = [r['R2'] for r in rf]
ax[0].semilogx(R2s, [r['n_front6_mean'] / 1e24 for r in rf], 'o-', color=C[0], label='front 6 sections')
ax[0].semilogx(R2s, [r['n_mean'] / 1e24 for r in rf], 's--', color=C[2], label='whole device')
ax[0].axhline(nz.net_transparency_density() / 1e24, color='k', ls=':', lw=0.8, label='net transparency, 1550 nm')
ax[0].set_xlabel('Rear-facet reflectivity $R_2$'); ax[0].set_ylabel(r'Mean $n$ (10$^{24}$ m$^{-3}$)')
ax[0].set_ylim(1.4, 1.75); ax[0].legend(fontsize=6, loc='upper right'); panel(ax[0], 'a')

ax[1].semilogx(R2s, [r['NF_A_dB'] for r in rf], 'o-', color=C[1])
ax[1].axhline(nz.QUANTUM_LIMIT_DB, color='k', ls='--', lw=0.8)
ax[1].set_xlabel('Rear-facet reflectivity $R_2$'); ax[1].set_ylabel('Chip-internal NF (dB)')
ax[1].set_ylim(0, 16); panel(ax[1], 'b')

g = res['ase_grid']
nb = [r['n_ase_bins'] for r in g]
ax[2].plot(nb, [r['NF_A_dB'] for r in g], 'o-', color=C[0])
ax[2].set_xscale('log'); ax[2].minorticks_off()
ax[2].set_xticks(nb); ax[2].set_xticklabels([str(b) for b in nb])
ax[2].axhline(nz.QUANTUM_LIMIT_DB, color='k', ls='--', lw=0.8)
ax[2].set_xlabel('Number of ASE wavelength bins'); ax[2].set_ylabel('Chip-internal NF (dB)')
ax[2].set_ylim(0, 16); panel(ax[2], 'c')

w = res['wavelength']
ax[3].plot([r['lam_nm'] for r in w], [r['NF_A_dB'] for r in w], 'o-', color=C[0])
ax[3].axhline(nz.QUANTUM_LIMIT_DB, color='k', ls='--', lw=0.8)
ax[3].set_xlabel('Signal wavelength (nm)'); ax[3].set_ylabel('Chip-internal NF (dB)')
ax[3].set_ylim(0, 20); panel(ax[3], 'd')
fig.tight_layout(); save(fig, 'FigV_validation', 'figures')

print(f'done in {time.time() - t0:.1f} s')
for r in res['rear_facet']:
    print(f"R2 {r['R2']:5.2f}  G_chip {r['G_chip_dB']:6.2f}  NF {r['NF_A_dB']:6.2f}  n_front {r['n_front6_mean']/1e24:.3f}"
          f"  netloss {r['sections_net_loss']}")
for r in res['ase_grid']:
    print(f"bins {r['n_ase_bins']:4d}  G_chip {r['G_chip_dB']:6.2f}  NF {r['NF_A_dB']:6.3f}  B {r['NF_B_dB']:6.3f}")
for r in res['wavelength']:
    print(f"lam {r['lam_nm']}  G_chip {r['G_chip_dB']:6.2f}  NF {r['NF_A_dB']:6.2f}  netloss {r['sections_net_loss']}  nsp_eff {r['nsp_eff_gain_weighted']:.2f}")
