import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COL1 = 3.5     # single-column width [in]
COL2 = 7.2     # double-column width [in]
C = ['#1b4f8a', '#c0392b', '#2e8b57', '#d68910', '#6c3483', '#566573', '#17a589']


def setup():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.labelsize': 8,
        'axes.titlesize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
        'ytick.labelsize': 7, 'axes.linewidth': 0.7, 'lines.linewidth': 1.2,
        'lines.markersize': 3.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True, 'axes.grid': True,
        'grid.color': '#dddddd', 'grid.linewidth': 0.5, 'legend.frameon': False,
        'savefig.dpi': 400, 'savefig.bbox': 'tight', 'figure.dpi': 150,
    })


def panel(ax, s):
    ax.text(-0.22, 1.03, s, transform=ax.transAxes, fontweight='bold', fontsize=9, va='bottom')


def save(fig, name, outdir):
    fig.savefig(f'{outdir}/{name}.png')
    fig.savefig(f'{outdir}/{name}.pdf')
    plt.close(fig)
