# Replica di "Zeta Spectral Triples" a 80 cifre (arXiv:2511.22755)
# Obiettivo: risolvere il cluster degenere degli autovalori minimi e
# chiudere la coda #45-50 della Tabella 6 (lam2 = 13, 14).
#
# Upgrade rispetto a replica_zst_hp.py (dps=35):
#  - dps = 85 (riportiamo a 80)
#  - Gauss-Legendre grado 60 per pannello di mezza oscillazione
#  - code oscillanti: integrazione per parti con derivate ESATTE
#    (Leibniz su w*R: polygamma per w, frazioni parziali per R)
#  - rifinitura zeri con findroot dopo bisezione

import json, time
from mpmath import (mp, mpf, mpc, sin, cos, sinh, sqrt, log, pi, psi,
                    quad, matrix, eigsy, legendre, binomial, factorial,
                    findroot, zetazero)
import numpy as np

mp.dps = 85
T0 = time.time()
def stamp(msg): print(f"[{time.time()-T0:7.0f}s] {msg}", flush=True)

# ---------- nodi GL raffinati ----------
def nodi_gl(g=60):
    x64, _ = np.polynomial.legendre.leggauss(g)
    nodi, pesi = [], []
    for x0 in x64:
        x = mpf(float(x0))
        for _ in range(6):
            Pn = legendre(g, x)
            dP = g * (x * Pn - legendre(g - 1, x)) / (x * x - 1)
            x = x - Pn / dP
        dP = g * (x * legendre(g, x) - legendre(g - 1, x)) / (x * x - 1)
        nodi.append(x); pesi.append(2 / ((1 - x * x) * dP * dP))
    return nodi, pesi
GLX, GLW = nodi_gl(60)

def von_mangoldt_dict(M):
    d = {}
    for p in range(2, M + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            pk, lp = p, log(p)
            while pk <= M:
                d[pk] = lp; pk *= p
    return d

K_IBP = 36
SEGNI = [-1, -1, +1, +1]; TRIG = [sin, cos, sin, cos]

def run(LAM2, N=120, S=mpf(2000)):
    L = log(LAM2); A = L / 2; b = 2 * A
    omega = [2 * pi * n / L for n in range(N + 1)]
    norm = [sqrt(1 / L)] + [sqrt(2 / L)] * N
    segno = [(-1) ** n for n in range(N + 1)]

    Av = [norm[n] * segno[n] * sinh(L / 4) / (mpf(1) / 4 + omega[n] ** 2)
          for n in range(N + 1)]

    # ---- primi ----
    stamp(f"lam2={LAM2}: termine dei primi")
    lamd = von_mangoldt_dict(LAM2)
    def corr(a, c, t):
        wa, wc = omega[a], omega[c]; tot = mpf(0)
        for Om in (wa + wc, wa - wc):
            if abs(Om) < mpf("1e-60"):
                tot += (2 * A - t) * cos(wa * t)
            else:
                tot += (sin(Om * (A - t) + wa * t) - sin(-Om * A + wa * t)) / Om
        return mpf("0.5") * norm[a] * norm[c] * tot
    PR = [[mpf(0)] * (N + 1) for _ in range(N + 1)]
    for k, lk in lamd.items():
        t = log(k); c0 = -2 * lk / sqrt(k)
        for a in range(N + 1):
            for c in range(a, N + 1):
                v = c0 * mpf("0.5") * (corr(a, c, t) + corr(c, a, t))
                PR[a][c] += v
                if a != c: PR[c][a] += v

    # ---- archimedeo: griglia principale ----
    stamp(f"lam2={LAM2}: griglia archimedea (digamma)")
    panw = pi / (2 * A); npan = int(S / panw) + 1
    nodi, base_int = [], []
    for kp in range(npan):
        lo, hi = kp * panw, min((kp + 1) * panw, S)
        c0, c1 = (hi + lo) / 2, (hi - lo) / 2
        for x, wgt in zip(GLX, GLW):
            s = c0 + c1 * x
            wv = psi(0, mpf("0.25") + mpf("0.5") * 1j * s).real - log(pi)
            nodi.append(s); base_int.append(wgt * c1 * sin(s * A) ** 2 * wv)
    stamp(f"  {len(nodi)} nodi pronti")

    # ---- derivate esatte di w in S (condivise) ----
    wder = [psi(0, mpf("0.25") + mpf("0.5") * 1j * S).real - log(pi)]
    for j in range(1, K_IBP + 1):
        wder.append(((mpc(0, "0.5")) ** j * psi(j, mpf("0.25") + mpf("0.5") * 1j * S)).real)

    def tail(om, quadrato):
        """coda di int_S^inf sin^2(sA) w(s) R(s) ds, R = 1/(s^2-om^2)^(1 o 2)."""
        # derivate esatte di R in S
        Rd = []
        if om == 0:
            for m in range(K_IBP + 1):
                e = (m + 2) if not quadrato else (m + 4)
                cst = factorial(m + 1) if not quadrato else factorial(m + 3) / 6
                Rd.append((-1) ** m * cst / S ** e)
        else:
            sm, sp = S - om, S + om
            for m in range(K_IBP + 1):
                f1 = (-1) ** m * factorial(m)
                if not quadrato:
                    Rd.append(f1 / (2 * om) * (sm ** (-(m + 1)) - sp ** (-(m + 1))))
                else:
                    f2 = (-1) ** m * factorial(m + 1)
                    Rd.append(f1 * (-1 / (4 * om ** 3)) * (sm ** (-(m + 1)) - sp ** (-(m + 1)))
                              + f2 / (4 * om ** 2) * (sm ** (-(m + 2)) + sp ** (-(m + 2))))
        # parte liscia: (1/2) int_S^inf w R ds  (u-sub, non oscillante)
        def rat_u(u):
            if not quadrato:
                return 1 / (1 / (u * u) - om * om) / (u * u)
            return 1 / (1 / (u * u) - om * om) ** 2 / (u * u)
        liscia = mpf("0.5") * quad(
            lambda u: (psi(0, mpf("0.25") + 1j / (2 * u)).real - log(pi)) * rat_u(u),
            [0, 1 / S])
        # parte oscillante: -(1/2) int_S^inf cos(bs) w R ds via IBP esatta
        I = mpf(0)
        for k in range(K_IBP + 1):
            phik = sum(binomial(k, j) * wder[j] * Rd[k - j] for j in range(k + 1))
            I += SEGNI[k % 4] * TRIG[k % 4](b * S) * phik / b ** (k + 1)
        return liscia - mpf("0.5") * I

    # ---- K e J in un passaggio ----
    stamp(f"lam2={LAM2}: somme K/J")
    Kv, Jv = {}, {}
    for n in range(N + 1):
        om2 = omega[n] ** 2
        sK = mpf(0); sJ = mpf(0)
        for i in range(len(nodi)):
            r = 1 / (nodi[i] ** 2 - om2)
            t = base_int[i] * r
            sK += t
            if n > 0: sJ += t * r
        Kv[n] = sK + tail(omega[n], False)
        if n > 0: Jv[n] = sJ + tail(omega[n], True)
        if n % 30 == 0: stamp(f"  K/J n={n}")

    AR = [[mpf(0)] * (N + 1) for _ in range(N + 1)]
    for a in range(N + 1):
        pa = omega[a] ** 2
        for c in range(a, N + 1):
            pb = omega[c] ** 2
            f = norm[a] * norm[c] * segno[a] * segno[c] / pi
            if a == c:
                v = f * (4 * Kv[a] + (4 * pa * Jv[a] if a > 0 else 0))
            else:
                v = f * 4 / (pa - pb) * (pa * Kv[a] - pb * Kv[c])
            AR[a][c] = v; AR[c][a] = v

    # ---- diagonalizzazione ----
    stamp(f"lam2={LAM2}: diagonalizzazione")
    M = matrix(N + 1)
    for a in range(N + 1):
        for c in range(N + 1):
            M[a, c] = 2 * Av[a] * Av[c] + PR[a][c] + AR[a][c]
    E, Q = eigsy(M)
    imin = min(range(N + 1), key=lambda i: E[i])
    ordinati = sorted(E[i] for i in range(N + 1))
    stamp("  5 autovalori minimi: " + ", ".join(mp.nstr(e, 4) for e in ordinati[:5]))
    xi = [Q[n, imin] for n in range(N + 1)]

    # ---- zeri di xi_hat ----
    def xhat(z):
        tot = mpf(0); sz = sin(z * A)
        for n in range(N + 1):
            if n == 0:
                tot += xi[0] * norm[0] * 2 * sz / z
            else:
                tot += xi[n] * norm[n] * segno[n] * 2 * z * sz / (z * z - omega[n] ** 2)
        return tot
    stamp(f"lam2={LAM2}: ricerca zeri")
    zeri = []
    zp, vp = mpf("0.5"), xhat(mpf("0.5"))
    z = zp + mpf("0.02")
    while z < 145:
        v = xhat(z)
        if vp * v < 0:
            lo, hi, vlo = zp, z, vp
            for _ in range(40):
                mid = (lo + hi) / 2; vm = xhat(mid)
                if vlo * vm < 0: hi = mid
                else: lo, vlo = mid, vm
            try:
                r = findroot(xhat, (lo + hi) / 2)
                zeri.append(r.real if hasattr(r, "real") else r)
            except Exception:
                zeri.append((lo + hi) / 2)
        zp, vp = z, v
        z += mpf("0.02")
    return zeri

# ---------- esecuzione ----------
print("zeri di Riemann di riferimento...")
GAMME = [zetazero(j).imag for j in range(1, 51)]
TAB6 = {13: {1: "2.44e-55", 10: "3.98e-37", 50: "2.04e-3"},
        14: {1: "1.07e-60", 10: "2.96e-42", 50: "4.78e-6"}}
risultati = {}
for LAM2 in (13, 14):
    zeri = run(LAM2)
    err = {}
    for j, g in enumerate(GAMME, start=1):
        err[j] = min(abs(zz - g) for zz in zeri)
    risultati[LAM2] = {j: float(e) for j, e in err.items()}
    print(f"\n=== lam2 = {LAM2} (dps 80): errori vs zeri di Riemann ===", flush=True)
    for j in (1, 5, 10, 20, 30, 40, 45, 48, 50):
        print(f"  zero #{j:>2}: nostro errore {float(err[j]):.3e}"
              f"   (paper: {TAB6[LAM2].get(j, '-')})", flush=True)

with open("replica_zst_80.json", "w") as f:
    json.dump(risultati, f, indent=1)
stamp("salvato: replica_zst_80.json")
