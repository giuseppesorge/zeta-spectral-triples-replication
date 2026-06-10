# Replica AD ALTA PRECISIONE di "Zeta Spectral Triples" (arXiv:2511.22755)
# mpmath dps=35, N=120, lambda^2 in {12, 13, 14} (Tabella 6 del paper).
#
# Bersaglio quantitativo: il pattern degli errori della Tabella 6 —
#   zero #50:  9.02e-2 (lam2=12),  2.04e-3 (lam2=13),  4.78e-6 (lam2=14)
# tutti sopra il pavimento numerico atteso (~1e-12), quindi verificabili.

import json, time
import numpy as np
from mpmath import (mp, mpf, sin, cos, sinh, sqrt, log, pi, digamma,
                    quad, diff, matrix, eigsy, legendre)

mp.dps = 35
T0 = time.time()

def nodi_gl(grado=30):
    """Nodi/pesi di Gauss-Legendre raffinati a precisione mpmath."""
    x64, w64 = np.polynomial.legendre.leggauss(grado)
    nodi, pesi = [], []
    for x0 in x64:
        x = mpf(float(x0))
        for _ in range(4):                       # Newton su P_n
            Pn = legendre(grado, x)
            dP = grado * (x * Pn - legendre(grado - 1, x)) / (x * x - 1)
            x = x - Pn / dP
        dP = grado * (x * legendre(grado, x) - legendre(grado - 1, x)) / (x * x - 1)
        nodi.append(x)
        pesi.append(2 / ((1 - x * x) * dP * dP))
    return nodi, pesi

GLX, GLW = nodi_gl(30)

def von_mangoldt_dict(M):
    d = {}
    for p in range(2, M + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            pk, lp = p, log(p)
            while pk <= M:
                d[pk] = lp
                pk *= p
    return d

def run(LAM2, N=120, S=mpf(1500)):
    L = log(LAM2); A = L / 2
    omega = [2 * pi * n / L for n in range(N + 1)]
    norm = [sqrt(1 / L)] + [sqrt(2 / L)] * N
    segno = [(-1) ** n for n in range(N + 1)]

    # ---------- polo (chiuso) ----------
    Av = [norm[n] * segno[n] * sinh(L / 4) / (mpf(1) / 4 + omega[n] ** 2)
          for n in range(N + 1)]

    # ---------- primi (chiuso) ----------
    lamd = von_mangoldt_dict(LAM2)
    def corr(a, b, t):
        wa, wb = omega[a], omega[b]
        tot = mpf(0)
        for Om in (wa + wb, wa - wb):
            if abs(Om) < mpf("1e-30"):
                tot += (2 * A - t) * cos(wa * t)
            else:
                tot += (sin(Om * (A - t) + wa * t) - sin(-Om * A + wa * t)) / Om
        return mpf("0.5") * norm[a] * norm[b] * tot
    print(f"[{time.time()-T0:6.0f}s] lam2={LAM2}: termine dei primi...")
    PR = [[mpf(0)] * (N + 1) for _ in range(N + 1)]
    for k, lk in lamd.items():
        t = log(k); c = -2 * lk / sqrt(k)
        for a in range(N + 1):
            for b in range(a, N + 1):
                v = c * mpf("0.5") * (corr(a, b, t) + corr(b, a, t))
                PR[a][b] += v
                if a != b: PR[b][a] += v

    # ---------- archimedeo: griglia globale + K(omega), J(omega) ----------
    print(f"[{time.time()-T0:6.0f}s] lam2={LAM2}: griglia archimedea...")
    panw = pi / (2 * A)                       # mezza oscillazione di sin^2
    npan = int(S / panw) + 1
    nodi, pesi, wval = [], [], []
    for kpan in range(npan):
        lo, hi = kpan * panw, min((kpan + 1) * panw, S)
        c0, c1 = (hi + lo) / 2, (hi - lo) / 2
        for x, wgt in zip(GLX, GLW):
            s = c0 + c1 * x
            nodi.append(s); pesi.append(wgt * c1)
            wval.append(digamma(mpf("0.25") + mpf("0.5") * 1j * s).real - log(pi))
    print(f"[{time.time()-T0:6.0f}s]   {len(nodi)} nodi; calcolo K, J...")
    sin2 = [sin(s * A) ** 2 for s in nodi]
    base_int = [pesi[i] * sin2[i] * wval[i] for i in range(len(nodi))]

    def code(rat, drat_ordini=6):
        """coda oltre S: parte liscia (u-sub) + oscillante (IBP)."""
        liscia = mpf("0.5") * quad(
            lambda u: (digamma(mpf("0.25") + 1j / (2 * u)).real - log(pi))
                      * rat(1 / u) / (u * u), [0, 1 / S])
        phi = lambda s: (digamma(mpf("0.25") + mpf("0.5") * 1j * s).real
                         - log(pi)) * rat(s)
        b = 2 * A; I = mpf(0); segni = [-1, -1, +1, +1]
        trig = [sin, cos, sin, cos]
        for k in range(drat_ordini):
            dk = diff(phi, S, k) if k else phi(S)
            I += segni[k % 4] * trig[k % 4](b * S) * dk / b ** (k + 1)
        return liscia - mpf("0.5") * I

    Kv, Jv = {}, {}
    for n in range(N + 1):
        om = omega[n]; om2 = om * om
        Kv[n] = (sum(base_int[i] / (nodi[i] ** 2 - om2) for i in range(len(nodi)))
                 + code(lambda s, o2=om2: 1 / (s * s - o2)))
        if n > 0:
            Jv[n] = (sum(base_int[i] / (nodi[i] ** 2 - om2) ** 2
                         for i in range(len(nodi)))
                     + code(lambda s, o2=om2: 1 / (s * s - o2) ** 2))
        if n % 30 == 0:
            print(f"[{time.time()-T0:6.0f}s]   K/J fino a n={n}")

    AR = [[mpf(0)] * (N + 1) for _ in range(N + 1)]
    for a in range(N + 1):
        pa = omega[a] ** 2
        for b in range(a, N + 1):
            pb = omega[b] ** 2
            f = norm[a] * norm[b] * segno[a] * segno[b] / pi
            if a == b:
                v = f * (4 * Kv[a] + (4 * pa * Jv[a] if a > 0 else 0))
            else:
                v = f * 4 / (pa - pb) * (pa * Kv[a] - pb * Kv[b])
            AR[a][b] = v; AR[b][a] = v

    # ---------- forma di Weil, autovettore minimo ----------
    print(f"[{time.time()-T0:6.0f}s] lam2={LAM2}: diagonalizzazione {N+1}x{N+1}...")
    M = matrix(N + 1)
    for a in range(N + 1):
        for b in range(N + 1):
            M[a, b] = 2 * Av[a] * Av[b] + PR[a][b] + AR[a][b]
    E, Q = eigsy(M)
    imin = min(range(N + 1), key=lambda i: E[i])
    print(f"[{time.time()-T0:6.0f}s]   autovalori minimi: "
          + ", ".join(mp.nstr(E[i], 3) for i in range(3)))
    xi = [Q[n, imin] for n in range(N + 1)]

    # ---------- zeri di xi_hat ----------
    def xhat(z):
        tot = mpf(0); sz = sin(z * A)
        for n in range(N + 1):
            if n == 0:
                tot += xi[0] * norm[0] * 2 * sz / z
            else:
                den = z * z - omega[n] ** 2
                tot += xi[n] * norm[n] * segno[n] * 2 * z * sz / den
        return tot
    print(f"[{time.time()-T0:6.0f}s] lam2={LAM2}: ricerca zeri...")
    zeri = []
    zprev, vprev = mpf("0.5"), xhat(mpf("0.5"))
    step = mpf("0.02")
    z = zprev + step
    while z < 145:
        v = xhat(z)
        if vprev * v < 0:
            lo, hi, vlo = zprev, z, vprev
            for _ in range(140):
                mid = (lo + hi) / 2; vm = xhat(mid)
                if vlo * vm < 0: hi = mid
                else: lo, vlo = mid, vm
            zeri.append((lo + hi) / 2)
        zprev, vprev = z, v
        z += step
    return zeri

# ---------------- esecuzione e confronto con la Tabella 6 ----------------
from mpmath import zetazero
print("zeri di Riemann di riferimento (alta precisione)...")
GAMME = [zetazero(j).imag for j in range(1, 51)]
TAB6 = {12: {1: "3.41e-50", 10: "2.99e-32", 50: "9.02e-2"},
        13: {1: "2.44e-55", 10: "3.98e-37", 50: "2.04e-3"},
        14: {1: "1.07e-60", 10: "2.96e-42", 50: "4.78e-6"}}

risultati = {}
for LAM2 in (12, 13, 14):
    zeri = run(LAM2)
    err = {}
    for j, g in enumerate(GAMME, start=1):
        vicino = min(zeri, key=lambda zz: abs(zz - g))
        err[j] = abs(vicino - g)
    risultati[LAM2] = {j: float(e) for j, e in err.items()}
    print(f"\n=== lam2 = {LAM2}: errori vs zeri di Riemann ===")
    for j in (1, 5, 10, 20, 30, 40, 45, 48, 50):
        atteso = TAB6[LAM2].get(j, "-")
        print(f"  zero #{j:>2}: nostro errore {float(err[j]):.3e}"
              f"   (paper: {atteso})")

with open("replica_zst_hp.json", "w") as f:
    json.dump(risultati, f, indent=1)
print(f"\n[{time.time()-T0:6.0f}s] salvato: replica_zst_hp.json")
