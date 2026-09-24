import pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
cn = pd.read_csv("cn_distribution_voronoi.csv", dtype={"rate": str})
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
g = cn[(cn["T"]==300)&(cn.rate=="2e10")&(cn.especie=="all")]
ax[0].bar(g.CN, g.fraccion, width=1, edgecolor="k", color="#1f77b4")
ax[0].set(xlim=(8,20), xticks=range(8,21,2), xlabel="Coordination Number (CN)", ylabel="Frequency")
ax[0].grid(axis="y", alpha=.4); ax[0].set_axisbelow(True)
for sp, col in [("Zr","#4c72b0"),("Cu","#dd8452"),("Ag","#55a868")]:
    g = cn[(cn["T"]==300)&(cn.rate=="2e10")&(cn.especie==sp)]
    ax[1].bar(g.CN, g.fraccion, width=1, color=col, alpha=.7, edgecolor="k", linewidth=.6, label=sp)
ax[1].set(xlim=(8,21), xticks=range(8,22,2), xlabel="Coordination Number (CN)", ylabel="Probability Density (PDF)")
ax[1].grid(alpha=.3); ax[1].set_axisbelow(True); ax[1].legend(frameon=True, fontsize=8, loc="upper right")
fig.suptitle("Voronoi (OVITO default), 2e10 K/s, 300 K", fontsize=9, y=1.0)
fig.tight_layout(); fig.savefig("fig6_replica_2e10_300K.png", dpi=150)
