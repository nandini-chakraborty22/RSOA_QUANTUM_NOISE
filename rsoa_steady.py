"""
Steady-state segmented travelling-wave model of a reflective SOA (RSOA).

Structure: front facet (z=0, reflectivity R1, fibre-coupled, seed in / signal out)
           rear facet  (z=L, reflectivity R2, high-reflection coated).
Signal and ASE are treated as incoherent photon-rate waves (phase-averaged over
cavity resonances); the coherent gain ripple is evaluated afterwards from the
single-pass gain with the Fabry-Perot reflection-amplifier formula.
Carrier density obeys Connelly (2001) Eq. (48) in each section.
"""
import numpy as np
import material as m

LAM_S = 1550e-9


class RSOA:
    def __init__(self, L=600e-6, R1=1e-4, R2=0.9, I=None, J=None,
                 eta_dB=-3.0, Kseg=40, n_ase=48, lam_s=LAM_S, p=m.P):
        self.p = p
        self.L, self.R1, self.R2 = L, R1, R2
        self.K = Kseg
        self.dz = L / Kseg
        self.V = p['W'] * p['d'] * L
        if I is None:
            I = J * p['W'] * L          # J in A/m^2 over the stripe area
        self.I = I
        self.eta = 10 ** (eta_dB / 10)
        self.lam_s = lam_s
        self.nu_s = m.c / lam_s
        lam = np.linspace(1.42e-6, 1.68e-6, n_ase)
        self.nu_a = m.c / lam
        self.dnu = np.abs(np.gradient(self.nu_a))

    # -------------------------------------------------------------------
    def solve(self, Pin_W, n_init=None, tol=1e-4, maxit=6000, ase=True):
        p, K, dz = self.p, self.K, self.dz
        G = p['Gamma']
        A = p['W'] * p['d']
        Jvol = self.I / (m.q * self.V)
        if n_init is None:
            lo, hi = 1e22, 1e26                      # n where R(n) = J/(qd): no-photon limit
            for _ in range(100):
                mid = np.sqrt(lo * hi)
                if m.recombination(mid, p) > Jvol: hi = mid
                else: lo = mid
            n = np.full(K, min(mid, 1.3e24))
        else:
            n = n_init.copy()
        Qin = self.eta * (1 - self.R1) * Pin_W / (m.h * self.nu_s)  # photons/s into chip
        M = len(self.nu_a)
        Sf = np.zeros(K + 1); Sb = np.zeros(K + 1)
        Af = np.zeros((K + 1, M)); Ab = np.zeros((K + 1, M))
        w = 0.5
        for it in range(maxit):
            gs, _ = m.gain_components(self.nu_s, n, p)
            ga, ge = m.gain_components(self.nu_a[None, :], n[:, None], p)
            alpha = m.alpha_loss(n, p)
            net_s = G * gs - alpha                       # (K,)
            net_a = G * ga - alpha[:, None]              # (K,M)
            Ts = np.exp(np.clip(net_s * dz, -50, 3))
            Ta = np.exp(np.clip(net_a * dz, -50, 3))
            # spontaneous source per section integrated (photons/s), per polarisation
            src = G * ge * self.dnu[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                Sa = np.where(np.abs(net_a) > 1e-9, src * (Ta - 1) / net_a, src * dz)
            if not ase:
                Sa = Sa * 0
            # iterate boundary conditions (multiple passes)
            for _ in range(6):
                Sf[0] = Qin + self.R1 * Sb[0]
                Af[0] = self.R1 * Ab[0]
                for k in range(K):
                    Sf[k + 1] = Sf[k] * Ts[k]
                    Af[k + 1] = Af[k] * Ta[k] + Sa[k]
                Sb[K] = self.R2 * Sf[K]
                Ab[K] = self.R2 * Af[K]
                for k in range(K - 1, -1, -1):
                    Sb[k] = Sb[k + 1] * Ts[k]
                    Ab[k] = Ab[k + 1] * Ta[k] + Sa[k]
            Smid = 0.5 * (Sf[:-1] + Sf[1:] + Sb[:-1] + Sb[1:])
            Amid = 0.5 * (Af[:-1] + Af[1:] + Ab[:-1] + Ab[1:])
            stim = G * gs * Smid / A + 2 * np.sum(G * ga * Amid, axis=1) / A
            # Connelly (2001) adaptive-step carrier update: F = dn/dt with frozen photon rates
            F = Jvol - m.recombination(n, p) - stim
            sgn = np.sign(F)
            if it == 0:
                eps = np.full(K, 0.05)
                prev = sgn.copy()
            flip = sgn != prev
            eps = np.where(flip, eps / 2, np.minimum(eps * 1.05, 0.05))
            prev = sgn
            n_new = np.clip(n * (1 + eps * sgn), 1e22, 1e26)
            rel = np.max(np.abs(F) / (m.recombination(n, p) + stim + Jvol))
            n = n_new
            if rel < tol and it > 20:
                break
        self.n = n
        self.Ts, self.net_s = Ts, net_s
        self._Amid = Amid
        hv = m.h * self.nu_s
        Pout = self.eta * (1 - self.R1) * Sb[0] * hv
        G_ff = Pout / Pin_W if Pin_W > 0 else np.nan
        # ASE out of the front facet, both polarisations
        ase_spec = 2 * self.eta * (1 - self.R1) * Ab[0] * m.h * self.nu_a   # W per bin
        Gs1 = np.prod(Ts)                                # single pass internal gain
        self.res = dict(Pout=Pout, G_ff=G_ff, Gs=Gs1, n=n, it=it,
                        ase_bins=ase_spec, ase_total=ase_spec.sum(),
                        ase_psd=ase_spec / self.dnu,     # W/Hz (both pols)
                        Sf=Sf.copy(), Sb=Sb.copy())
        return self.res

    # -------------------------------------------------------------------
    def ripple_dB(self):
        """Peak-to-valley gain ripple of a reflection-type FP amplifier."""
        r = np.sqrt(self.R1 * self.R2) * self.res['Gs']
        if r >= 1:
            return np.inf
        Gs, R1, R2 = self.res['Gs'], self.R1, self.R2
        def Gr(s2):
            num = (np.sqrt(R1) - np.sqrt(R2) * Gs) ** 2 + 4 * np.sqrt(R1 * R2) * Gs * s2
            den = (1 - np.sqrt(R1 * R2) * Gs) ** 2 + 4 * np.sqrt(R1 * R2) * Gs * s2
            return num / den
        return 10 * np.log10(Gr(0.0) / Gr(1.0))

    def lasing_margin(self):
        """sqrt(R1 R2) * Gs ; >= 1 means the cavity reaches lasing threshold."""
        return np.sqrt(self.R1 * self.R2) * self.res['Gs']

    def tau_eff(self):
        """Effective (differential) carrier lifetime averaged along the device [s]."""
        p, n = self.p, self.n
        A = p['W'] * p['d']
        Sf, Sb = self.res['Sf'], self.res['Sb']
        Smid = 0.5 * (Sf[:-1] + Sf[1:] + Sb[:-1] + Sb[1:])
        dn_ = 1e18
        g1, _ = m.gain_components(self.nu_s, n, p)
        g2, _ = m.gain_components(self.nu_s, n + dn_, p)
        a_diff = (g2 - g1) / dn_
        rate = m.dR_dn(n, p) + p['Gamma'] * a_diff * Smid / A
        return 1.0 / np.mean(rate)

    def tau_diff_ase(self):
        """Differential carrier lifetime including ASE-stimulated recombination but
        excluding the injected signal (the signal term is handled by the dynamic model)."""
        p, n = self.p, self.n
        A = p['W'] * p['d']
        dn_ = 1e18
        g1, _ = m.gain_components(self.nu_a[None, :], n[:, None], p)
        g2, _ = m.gain_components(self.nu_a[None, :], (n + dn_)[:, None], p)
        rate = m.dR_dn(n, p) + 2 * np.sum(p['Gamma'] * (g2 - g1) / dn_ * self._Amid, axis=1) / A
        return 1.0 / np.mean(rate)
