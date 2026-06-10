# REPLICA INDIPENDENTE (a bassa precisione) di "Zeta Spectral Triples"
# (Connes-Consani-Moscovici, arXiv:2511.22755, novembre 2025).
#
# Costruzione (estratta dal paper):
#  - spazio: L^2([lambda^-1, lambda], du/u); in coordinate y = log u,
#    l'intervallo [-L/2, L/2] con L = 2 log(lambda);
#  - sul sottospazio PARI (funzioni invarianti per u -> 1/u) con base
#    coseno c_n(y), si costruisce la FORMA QUADRATICA DI WEIL troncata:
#        QW = polo - primi + archimedeo
#    (la stessa del nostro esperimento positivita_weil.py: entrano solo
#    i primi k <= lambda^2 perche' il supporto delle correlazioni e' [-L,L]);
#  - si prende l'autovettore xi di autovalore MINIMO epsilon_N (la funzione
#    test "piu' vicina a violare la positivita' di Weil");
#  - TESI: la trasformata di Fourier xi_hat(z) (funzione intera) ha zeri
#    reali che convergono agli zeri non banali di zeta. Dimostrare la
#    convergenza per N, lambda -> infinito equivarrebbe a RH.
#
# Il paper usa N=120 e 200 cifre di precisione (errori fino a 1e-60).
# Qui: N piu' piccolo, float64 + quadratura -> replica STRUTTURALE:
# verifichiamo che gli zeri emergano e a quante cifre coincidono.

import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpmath import zetazero

# ---------------- parametri ----------------
LAM2 = 13          # lambda^2: primi inclusi k <= 13 (come nel paper, sqrt(13))
N    = 60          # modi coseno 0..N
L    = np.log(LAM2)  # L = 2 log lambda = log(lambda^2)
A    = L / 2.0
omega = 2 * np.pi * np.arange(N + 1) / L
norm  = np.where(np.arange(N + 1) == 0, np.sqrt(1 / L), np.sqrt(2 / L))

# ---------------- termine del polo ----------------
# A_pm(n) = integrale di c_n(y) e^{\pm y/2} dy = norm*(-1)^n*sinh(L/4)/(1/4+omega^2)
Avec = norm * ((-1.0) ** np.arange(N + 1)) * np.sinh(L / 4) / (0.25 + omega ** 2)
POLO = 2.0 * np.outer(Avec, Avec)

# ---------------- termine dei primi ----------------
def von_mangoldt(M):
    lam = np.zeros(M + 1)
    for p in range(2, M + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            pk = p
            while pk <= M:
                lam[pk] = np.log(p)
                pk *= p
    return lam

lam = von_mangoldt(LAM2)

def correlazione(t):
    """q_ab(t) = int_{-A}^{A-t} c_a(y+t) c_b(y) dy, in forma chiusa esatta."""
    Q = np.zeros((N + 1, N + 1))
    for a_ in range(N + 1):
        wa = omega[a_]
        for b_ in range(N + 1):
            wb = omega[b_]
            tot = 0.0
            for Om in (wa + wb, wa - wb):
                if abs(Om) < 1e-14:
                    tot += (2 * A - t) * np.cos(wa * t)
                else:
                    tot += (np.sin(Om * (A - t) + wa * t)
                            - np.sin(-Om * A + wa * t)) / Om
            Q[a_, b_] = 0.5 * norm[a_] * norm[b_] * tot
    return Q

PRIMI = np.zeros((N + 1, N + 1))
contributi_primi = {}                       # per la decomposizione per primo
for k in range(2, LAM2 + 1):
    if lam[k] > 0:
        Q = correlazione(np.log(k))
        contrib = -2.0 * lam[k] / np.sqrt(k) * 0.5 * (Q + Q.T)
        PRIMI += contrib
        contributi_primi[k] = contrib

# ---------------- termine archimedeo ----------------
# Arch(a,b) = (1/pi) int_0^inf chat_a(s) chat_b(s) w(s) ds,
# w(s) = Re psi(1/4 + i s/2) - log(pi);  chat_n in forma chiusa.
def chat(n, s):
    s = np.asarray(s, dtype=float)
    if n == 0:
        out = np.sqrt(1 / L) * 2 * np.sin(s * A) / s
        return np.where(np.abs(s) < 1e-12, np.sqrt(1 / L) * L, out)
    num = 2 * s * np.sin(s * A) * (-1.0) ** n
    den = s ** 2 - omega[n] ** 2
    out = np.sqrt(2 / L) * num / den
    # limite rimovibile in s = omega_n: chat = norm * A * (-1)^n * ... -> sqrt(2/L)*A...
    sing = np.abs(den) < 1e-9
    return np.where(sing, np.sqrt(2 / L) * A, out)

# griglia di quadratura: densa vicino alle risonanze, coda fino a S
S = 4000.0
s1 = np.linspace(1e-6, omega[-1] * 2, 120_000)
s2 = np.geomspace(omega[-1] * 2, S, 60_000)
sgrid = np.concatenate([s1, s2[1:]])
w = digamma(0.25 + 0.5j * sgrid).real - np.log(np.pi)
C = np.vstack([chat(n, sgrid) for n in range(N + 1)])
ARCH = (1 / np.pi) * np.einsum("ak,bk,k->ab", C, C * w[None, :],
                                np.gradient(sgrid))
# correzione di coda analitica per s > S (sin^2 -> 1/2, w(s) ~ log(s/2pi)):
segni = (-1.0) ** np.arange(N + 1)
coda = (1 / np.pi) * 2 * (np.log(S / (2 * np.pi)) + 1) / S
ARCH += coda * np.outer(norm * segni, norm * segni)

# ---------------- forma di Weil e autovettore minimo ----------------
QW = POLO + PRIMI + ARCH
QW = 0.5 * (QW + QW.T)
evals, evecs = eigh(QW)
print("autovalori minimi della forma di Weil troncata:")
print("  ", np.array2string(evals[:6], precision=3))
print("  (positivita' di Weil: tutti >= 0 a meno del rumore numerico)")
xi = evecs[:, 0]

# ---------------- zeri di xi_hat vs zeri di Riemann ----------------
def xi_hat(z):
    return sum(xi[n] * chat(n, z) for n in range(N + 1))

zs = np.linspace(0.5, 42, 60_000)
vals = xi_hat(zs)
zeri = []
for i in range(len(zs) - 1):
    if vals[i] * vals[i + 1] < 0:
        lo, hi = zs[i], zs[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if xi_hat(np.array([mid]))[0] * xi_hat(np.array([lo]))[0] < 0:
                hi = mid
            else:
                lo = mid
        zeri.append(0.5 * (lo + hi))

print("\nzeri di Riemann vs zero piu' vicino di xi_hat (replica):")
print("(gli altri zeri di xi_hat sono i modi del reticolo 2pi*k/L dell'intervallo,")
print(" presenti anche nello spettro dell'operatore del paper)")
gammas = [float(zetazero(j).imag) for j in range(1, 12)]
zeri_arr = np.array(zeri)
for g in gammas:
    if g < zeri_arr.max():
        z = zeri_arr[np.argmin(np.abs(zeri_arr - g))]
        print(f"  Riemann: {g:>12.6f}   replica: {z:>12.6f}   errore: {abs(z-g):.2e}")

# ---------------- grafici ----------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
ax = axes[0]
ax.plot(zs, vals / np.max(np.abs(vals)), lw=0.9, color="#1a1a72",
        label="ξ̂(z) (replica, normalizzata)")
for g in gammas:
    ax.axvline(g, color="crimson", lw=0.8, ls="--", alpha=0.7)
ax.axhline(0, color="gray", lw=0.6)
ax.set_title(f"Replica di Zeta Spectral Triples (λ²={LAM2}, N={N}, float64)\n"
             "linee rosse = zeri veri di ζ; ξ̂ costruita SOLO con i primi ≤ 13")
ax.set_xlabel("z"); ax.legend(fontsize=9)

ax = axes[1]
ks = sorted(contributi_primi)
vals_k = [float(xi @ contributi_primi[k] @ xi) for k in ks]
polo_v = float(xi @ POLO @ xi); arch_v = float(xi @ ARCH @ xi)
labels = ["polo"] + [str(k) for k in ks] + ["arch"]
heights = [polo_v] + vals_k + [arch_v]
colors = ["#2e7d32"] + ["#b71c1c"] * len(ks) + ["#e65100"]
ax.bar(labels, heights, color=colors)
ax.axhline(0, color="gray", lw=0.6)
ax.set_title("La battaglia sulla funzione minimizzante ξ:\n"
             "contributo del polo, di ogni potenza di primo, e del posto ∞")
ax.set_ylabel("contributo a W(ξ)")
plt.tight_layout()
plt.savefig("replica_zst.png", dpi=120)
print("\nsalvato: replica_zst.png")
print(f"W(xi) totale = {float(xi @ QW @ xi):.3e}  (= autovalore minimo)")
