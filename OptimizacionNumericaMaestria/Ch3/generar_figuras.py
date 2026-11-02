#!/usr/bin/env python3
"""Figuras de la presentación del Capítulo 3 (Brunton, Optimization Bootcamp).

Genera en images/:
  fig_3_1_politopo.png    -- región factible, conjuntos de nivel y vértices (Figs. 3.1/3.3)
  fig_3_5_pivoteo.png     -- camino simplex de vértice en vértice (Fig. 3.5)
  fig_3_8_relajacion.png  -- relajación LP frente a la solución entera (Panel 3.3)

Usa el backend Agg para correr sin ventana en el servidor.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")
os.makedirs(OUT, exist_ok=True)

UNAH = "#00529B"
ACENTO = "#00876E"
GRIS = "#9AA0A6"


# ---------------------------------------------------------------------------
# Figura 1: región factible del problema de los microscopios + conjuntos de nivel
# ---------------------------------------------------------------------------
def figura_politopo():
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)

    # región factible: (0,0), (50,0), (20,30), (0,30)
    poly = np.array([[0, 0], [50, 0], [20, 30], [0, 30]])
    ax.add_patch(plt.Polygon(poly, closed=True, facecolor=ACENTO, alpha=0.16,
                              edgecolor=ACENTO, linewidth=2.0, zorder=2))

    # conjuntos de nivel de f(x) = 100 x1 + 400 x2
    x1 = np.linspace(-4, 56, 400)
    for k in range(0, 14001, 2000):
        x2 = (k - 100 * x1) / 400.0
        ax.plot(x1, x2, color=GRIS, linewidth=0.9, alpha=0.85, zorder=1)

    # vértices
    verts = [((0, 0), "0"), ((50, 0), "5000"),
             ((20, 30), "14000"), ((0, 30), "12000")]
    for (vx, vy), val in verts:
        ax.plot([vx], [vy], "o", color=UNAH, markersize=8, zorder=4)
        ax.annotate(f"({vx}, {vy})\nf = {val}", (vx, vy),
                    textcoords="offset points", xytext=(9, 8),
                    fontsize=9, color=UNAH, weight="bold")

    ax.plot([0, 0], [0, 30], color="#444444", linewidth=1.4, linestyle="--", zorder=3)
    ax.plot([50, 20], [0, 30], color="#444444", linewidth=1.4, linestyle="--", zorder=3)
    ax.plot([0, 20], [30, 30], color="#444444", linewidth=1.4, linestyle="--", zorder=3)

    ax.set_xlim(-6, 60)
    ax.set_ylim(-6, 40)
    ax.set_xlabel("$x_1$  (microscopios de enseñanza)", fontsize=10)
    ax.set_ylabel("$x_2$  (microscopios de investigación)", fontsize=10)
    ax.set_title("Región factible (polítopo convexo) y conjuntos de nivel de $f$",
                 fontsize=11, color=UNAH)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.grid(alpha=0.18, linestyle=":")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_3_1_politopo.png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figura 2: el camino del simplex, de vértice en vértice
# ---------------------------------------------------------------------------
def figura_pivoteo():
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)

    poly = np.array([[0, 0], [50, 0], [20, 30], [0, 30]])
    ax.add_patch(plt.Polygon(poly, closed=True, facecolor=ACENTO, alpha=0.13,
                              edgecolor=ACENTO, linewidth=2.0, zorder=2))

    camino = [((0, 0), (0, 30)), ((0, 30), (20, 30))]
    for i, (p0, p1) in enumerate(camino, start=1):
        ax.annotate("", xy=p1, xytext=p0,
                    arrowprops=dict(arrowstyle="-|>", color=UNAH, linewidth=2.6,
                                    mutation_scale=20, shrinkA=9, shrinkB=9), zorder=5)
        mid = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        ax.text(mid[0] + (2 if i == 2 else -7), mid[1] + (3.0 if i == 2 else 2.5),
                f"pivote {i}", color=UNAH, fontsize=10, weight="bold", zorder=6)

    etiquetas = [((0, 0), "vértice inicial\n$f = 0$"),
                 ((0, 30), "$f = 12000$"),
                 ((20, 30), "óptimo\n$f = 14000$")]
    for (vx, vy), txt in etiquetas:
        ax.plot([vx], [vy], "o", color=UNAH, markersize=9, zorder=7)
        off = (12, -26) if (vx, vy) == (0, 0) else ((-74, 8) if (vx, vy) == (0, 30) else (10, 4))
        ax.annotate(txt, (vx, vy), textcoords="offset points", xytext=off,
                    fontsize=9.5, color=UNAH, weight="bold", zorder=7)

    ax.set_xlim(-8, 58)
    ax.set_ylim(-8, 40)
    ax.set_xlabel("$x_1$", fontsize=11)
    ax.set_ylabel("$x_2$", fontsize=11)
    ax.set_title("Simplex de Dantzig: se camina de vértice en vértice", fontsize=11, color=UNAH)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.grid(alpha=0.18, linestyle=":")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_3_5_pivoteo.png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figura 3: relajación LP frente a la solución entera (gap de relajamiento)
# ---------------------------------------------------------------------------
def figura_relajacion():
    # knapsack binario de 2 ítems: 3x1 + 5x2 <= 6,  max 4x1 + 7x2
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=200)

    xs = np.linspace(0, 1, 400)
    ax.fill_between(xs, 0, np.minimum(1.0, (6 - 3 * xs) / 5.0),
                    color=ACENTO, alpha=0.15, edgecolor="none", zorder=1)
    ax.plot(xs, np.minimum(1.0, (6 - 3 * xs) / 5.0), color=ACENTO, linewidth=2.0, zorder=3)

    # puntos enteros factibles
    enteros = [(0, 0), (1, 0), (0, 1)]
    for (a, b) in enteros:
        ax.plot([a], [b], "s", color="#444444", markersize=9, zorder=6)
    ax.plot([0], [1], "s", color="#444444", markersize=9, zorder=6)
    ax.annotate("enteros factibles\n$[(0,0),(1,0),(0,1)]$", (0.06, 0.55),
                fontsize=9, color="#444444")

    # óptimo LP (1/3, 1) con valor 8.333
    ax.plot([1 / 3], [1], "*", color="#C62828", markersize=18, zorder=7)
    ax.annotate("óptimo de la relajación LP\n$(1/3,\\,1)$,  $f = 8.33$", (1 / 3, 1),
                textcoords="offset points", xytext=(14, -34), fontsize=9.5,
                color="#C62828", weight="bold",
                arrowprops=dict(arrowstyle="->", color="#C62828", linewidth=1.2))

    # óptimo entero (0,1) con valor 7
    ax.plot([0], [1], "s", color="#1B5E20", markersize=13, zorder=8)
    ax.annotate("óptimo entero\n$(0,1)$,  $f = 7$", (0, 1),
                textcoords="offset points", xytext=(-6, 26), fontsize=9.5,
                color="#1B5E20", weight="bold", ha="center")

    ax.set_xlim(-0.12, 1.15)
    ax.set_ylim(-0.12, 1.30)
    ax.set_xlabel("$x_1$", fontsize=11)
    ax.set_ylabel("$x_2$", fontsize=11)
    ax.set_title("Relajación LP vs. óptimo entero: gap de relajamiento",
                 fontsize=11, color=UNAH)
    ax.grid(alpha=0.2, linestyle=":")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_3_8_relajacion.png"))
    plt.close(fig)


if __name__ == "__main__":
    figura_politopo()
    figura_pivoteo()
    figura_relajacion()
    for f in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, f)
        print(f"{f:<30} {os.path.getsize(p):>8} bytes")
