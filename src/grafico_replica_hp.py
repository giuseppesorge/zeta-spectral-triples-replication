# Figura riassuntiva della replica ad alta precisione:
# errore per indice di zero, per i tre valori di lambda^2, con i tre
# punti dichiarati nella Tabella 6 del paper sovrapposti.
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# i JSON vengono cercati prima nella directory corrente (run appena
# rigenerati), poi in results/ accanto a questo script (dati pubblicati)
RESULTS = Path(__file__).resolve().parent.parent / "results"
def carica(nome):
    p = Path(nome) if Path(nome).exists() else RESULTS / nome
    return json.load(open(p))

ris35 = carica("replica_zst_hp.json")
ris80 = carica("replica_zst_80.json")
paper = {"12": {1: 3.41e-50, 10: 2.99e-32, 50: 9.02e-2},
         "13": {1: 2.44e-55, 10: 3.98e-37, 50: 2.04e-3},
         "14": {1: 1.07e-60, 10: 2.96e-42, 50: 4.78e-6}}
colori = {"12": "#1a6b1a", "13": "#1a1a72", "14": "#b71c1c"}

fig, ax = plt.subplots(figsize=(11, 7))
for lam2, errs in ris35.items():
    idx = sorted(int(k) for k in errs)
    ax.semilogy(idx, [errs[str(j)] for j in idx], "--", lw=0.8, alpha=0.5,
                color=colori[lam2], label=f"replica 35 cifre, λ² = {lam2}")
for lam2, errs in ris80.items():
    idx = sorted(int(k) for k in errs)
    ax.semilogy(idx, [errs[str(j)] for j in idx], "o-", ms=3, lw=1.4,
                color=colori[lam2], label=f"replica 80 cifre, λ² = {lam2}")
for lam2 in ris80:
    pj = sorted(paper[lam2])
    ax.semilogy(pj, [paper[lam2][j] for j in pj], "*", ms=16,
                color=colori[lam2], markeredgecolor="black",
                label=f"paper (Tab. 6, 200 cifre) λ² = {lam2}")
ax.set_xlabel("indice dello zero di Riemann")
ax.set_ylabel("errore |zero replicato − zero vero|")
ax.set_title("Replica indipendente di Zeta Spectral Triples (arXiv:2511.22755)\n"
             "N = 120; la legge di scala in λ del paper è confermata; "
             "sotto il pavimento i valori del paper richiedono 200 cifre")
ax.legend(fontsize=9)
ax.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig("replica_zst_hp.png", dpi=120)
print("salvato: replica_zst_hp.png")
