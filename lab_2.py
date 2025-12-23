import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# ВХІДНІ ДАНІ
# Варіант №9: lam = 3.3 см, ap x bp = 12x12 см
# -----------------------------
lam = 0.033   # довжина хвилі, м (3.3 см)
ap  = 0.12    # розмір розкриву в H-площині, м
bp  = 0.12    # розмір розкриву в E-площині, м

NPTS = 5000
theta_deg = np.linspace(0, 90, NPTS)
theta = np.deg2rad(theta_deg)

LEVEL = 0.707  # рівень половинної потужності

# -----------------------------
# Допоміжні функції
# -----------------------------
def sinc_safe(x: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """sinc(x) = sin(x)/x з обробкою околу нуля."""
    out = np.ones_like(x)
    m = np.abs(x) > eps
    out[m] = np.sin(x[m]) / x[m]
    return out

def hpbw(theta_deg: np.ndarray, F_norm: np.ndarray, level: float = LEVEL):
    """ШГП на рівні level: повертає (2*θ_edge, θ_edge)."""
    hit = np.flatnonzero(F_norm <= level)
    if hit.size == 0:
        return None, None
    edge = float(theta_deg[hit[0]])
    return 2 * edge, edge

def local_peaks(theta_deg: np.ndarray, F_norm: np.ndarray, start_deg: float = 0.0):
    """Локальні максимуми F_norm починаючи з кута start_deg."""
    start_i = int(np.searchsorted(theta_deg, start_deg, side="left"))
    th = []
    val = []
    for i in range(max(start_i, 1), len(F_norm) - 1):
        if F_norm[i] > F_norm[i - 1] and F_norm[i] > F_norm[i + 1]:
            th.append(theta_deg[i])
            val.append(F_norm[i])
    return np.array(th), np.array(val)

def zeros_E(lam: float, bp: float, nmax: int = 10):
    """sin(xE)=0 => θ = arcsin(nλ/bp)."""
    z = []
    for n in range(1, nmax):
        v = n * lam / bp
        if v < 1:
            z.append(np.rad2deg(np.arcsin(v)))
    return np.array(z)

def zeros_H(lam: float, ap: float, nmax: int = 10):
    """cos(xH)=0 => xH = π/2(2n+1) => θ = arcsin((2n+1)λ/(2ap))."""
    z = []
    for n in range(1, nmax):
        v = (2 * n + 1) * lam / (2 * ap)
        if v < 1:
            z.append(np.rad2deg(np.arcsin(v)))
    return np.array(z)

def plot_plane(title: str, theta_deg: np.ndarray, base: np.ndarray, elem: np.ndarray,
               total_norm: np.ndarray, zeros: np.ndarray, edge: float, HPBW: float,
               sl_th: np.ndarray, sl_val: np.ndarray, base_lbl: str):
    plt.figure(figsize=(10, 6))
    plt.plot(theta_deg, base, label=base_lbl, linewidth=1.3)
    plt.plot(theta_deg, elem, label="F1(θ) = (1+cosθ)/2", linewidth=1.3)
    plt.plot(theta_deg, total_norm, label="F(θ) — повна ДС (нормована)", linewidth=1.5)

    plt.axhline(LEVEL, color="gray", linestyle="--", linewidth=0.8)
    if edge is not None and HPBW is not None:
        plt.scatter(edge, LEVEL, color="red",
                    label=f"ШГП ≈ {HPBW:.2f}° (права межа {edge:.2f}°)")

    if zeros.size:
        plt.plot(zeros, np.zeros_like(zeros), "ko", label="Нулі")
        for z in zeros:
            plt.scatter(z, 0, color="black")

    if sl_th.size:
        for a, v in zip(sl_th, sl_val):
            plt.scatter(a, v, color="magenta")
        sll1_db = 20 * np.log10(sl_val[0])
        plt.scatter(sl_th[0], sl_val[0], color="cyan",
                    label=f"1-й бічний пелюсток ≈ {sll1_db:.1f} дБ")

    plt.title(title)
    plt.xlabel("Кут θ, градуси")
    plt.ylabel("|F(θ)| (нормована)")
    plt.xlim(0, 90)
    plt.ylim(0, 1.05)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

# -----------------------------------------------------------
# Нулі за методичкою
# -----------------------------------------------------------
zE = zeros_E(lam, bp)
zH = zeros_H(lam, ap)

# ============================================================
#                     ПЛОЩИНА E
# ============================================================
xE = np.pi * bp / lam * np.sin(theta)
F_CE = sinc_safe(xE)
F1   = (1 + np.cos(theta)) / 2
F_E  = F_CE * F1
F_E_norm = np.abs(F_E) / np.max(np.abs(F_E))

HPBW_E, edge_E = hpbw(theta_deg, F_E_norm)
th_sl_E, val_sl_E = local_peaks(theta_deg, F_E_norm, start_deg=float(zE[0]) if zE.size else 0.0)

plot_plane(
    title="ДС пірамідального рупору в площині E (варіант 9)",
    theta_deg=theta_deg,
    base=F_CE,
    elem=F1,
    total_norm=F_E_norm,
    zeros=zE,
    edge=edge_E,
    HPBW=HPBW_E,
    sl_th=th_sl_E,
    sl_val=val_sl_E,
    base_lbl="F_CE(θ) — основний множник"
)

# ============================================================
#                     ПЛОЩИНА H
# ============================================================
xH = np.pi * ap / lam * np.sin(theta)
den = 1 - (2 * ap / lam * np.sin(theta))**2
den[np.isclose(den, 0)] = 1e-9

F_CH = np.cos(xH) / den
F_H  = F_CH * F1
F_H_norm = np.abs(F_H) / np.max(np.abs(F_H))

HPBW_H, edge_H = hpbw(theta_deg, F_H_norm)
th_sl_H, val_sl_H = local_peaks(theta_deg, F_H_norm, start_deg=float(zH[0]) if zH.size else 0.0)

plot_plane(
    title="ДС пірамідального рупору в площині H (варіант 9)",
    theta_deg=theta_deg,
    base=F_CH,
    elem=F1,
    total_norm=F_H_norm,
    zeros=zH,
    edge=edge_H,
    HPBW=HPBW_H,
    sl_th=th_sl_H,
    sl_val=val_sl_H,
    base_lbl="F_CH(θ) — основний множник"
)

# ============================================================
#                ТАБЛИЦІ ДЛЯ ЛАБОРАТОРНОЇ
# ============================================================
print("\nТабл. 1 – Значення нульових кутів (θmin)")
print("-------------------------------------------------")
print("| № | θmin_H | θmin_E |")
for i in range(max(len(zH), len(zE))):
    znH = f"{zH[i]:.2f}" if i < len(zH) else "  -  "
    znE = f"{zE[i]:.2f}" if i < len(zE) else "  -  "
    print(f"| {i+1} | {znH:7} | {znE:7} |")
print("-------------------------------------------------")

print("Табл. 2 – Значення максимальних кутів E (бічні пелюстки)")
print("-------------------------------------------------")
print("| № | θmax_E | FE(θ)  |")
for i in range(len(th_sl_E)):
    print(f"| {i+1} | {th_sl_E[i]:7.2f} | {val_sl_E[i]:7.3f} |")
print("-------------------------------------------------")

print("Табл. 3 – Значення максимальних кутів H (бічні пелюстки)")
print("-------------------------------------------------")
print("| № | θmax_H | FH(θ)  |")
for i in range(len(th_sl_H)):
    print(f"| {i+1} | {th_sl_H[i]:7.2f} | {val_sl_H[i]:7.3f} |")
print("-------------------------------------------------")

plt.show()