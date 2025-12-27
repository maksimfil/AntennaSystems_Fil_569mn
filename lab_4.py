import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# (Варіант 9)
# lambda_cm = 4.1 cm
# h = 14 cm
# ============================================================

HP_LEVEL = 1 / np.sqrt(2)
EPS = 1e-12

# ------------ 1) Вхідні дані -------------
eps_r = 2.2
tg_delta = 2e-4
l_cm      = 23.7

d_max_cm = 2.9
d_min_cm = 1.7
d_cp_cm  = 2.3

lambda_cm = 4.1
h_cm      = 14.0

lam = lambda_cm / 100
h   = h_cm / 100
l   = l_cm / 100

xi = 1 + lam / (2 * l)
print("xi =", xi)

# ------------ 2) Кутова сітка --------------
theta_deg = np.linspace(0, 90, 2001)
theta = np.deg2rad(theta_deg)

# ------------ 3) Функції ДС -----------------
def norm(x):
    x = np.abs(x)
    m = np.max(x)
    return x / m if m > 0 else x

def F_B(theta):
    k = np.pi * (l / lam)
    X = xi - np.cos(theta)

    num = np.sin(k * X)
    den = X * np.sin(k * (xi - 1))

    out = np.zeros_like(theta, dtype=float)
    mask = np.abs(den) > EPS
    out[mask] = num[mask] / den[mask]
    out[~mask] = 1.0
    return out

def F_C(theta):
    return np.cos(np.pi * h / lam * np.sin(theta))

def find_zeros(f):
    return list(np.where(f[:-1] * f[1:] < 0)[0] + 1)

def find_maxima(f):
    a = np.abs(f)
    return list(np.where((a[1:-1] > a[:-2]) & (a[1:-1] > a[2:]))[0] + 1)

def level_cross_first_angle(f, level=HP_LEVEL):
    fn = norm(f)
    idxs = np.where((fn[:-1] >= level) & (fn[1:] <= level))[0]
    if idxs.size == 0:
        return None

    i = int(idxs[0] + 1)
    x1, x2 = theta_deg[i - 1], theta_deg[i]
    y1, y2 = fn[i - 1], fn[i]
    if abs(y2 - y1) < EPS:
        return float(x2)

    x = x1 + (level - y1) * (x2 - x1) / (y2 - y1)
    return float(x)

def print_table(title, idxs, values=None, n=10):
    print("\n" + title)
    print("--------------------------------------")
    print("| № |   θ (град)   |   Значення      |")
    print("--------------------------------------")
    for j, idx in enumerate(idxs[:n], 1):
        ang = float(theta_deg[idx])
        if values is None:
            print(f"| {j:2d} |   {ang:.4f}    |     0.0000      |")
        else:
            print(f"| {j:2d} |   {ang:.4f}    |    {float(values[idx]):.4f}    |")
    print("--------------------------------------")

def annotate_zeros(ax, idxs, prefix, y_shift, color, marker, text_color, n=10):
    for i, idx in enumerate(idxs[:n], 1):
        x = theta_deg[idx]
        ax.scatter(x, 0, color=color, s=45, marker=marker)
        ax.text(x, y_shift, f"{prefix}{i}", ha="center",
                transform=ax.get_xaxis_transform(), color=text_color)

def annotate_maxima(ax, idxs, fn, prefix, color, marker, n=10, dx=1.0, dy=0.03):
    for i, idx in enumerate(idxs[:n], 1):
        x = theta_deg[idx]
        y = fn[idx]
        ax.scatter(x, y, color=color, marker=marker, s=70)
        ax.text(x + dx, y + dy, f"{prefix}{i}", color=color)

def annotate_level(ax, x, level, color, dy=0.03):
    if x is None:
        return
    ax.scatter(x, level, color=color, s=70, zorder=5)
    ax.text(x + 1, level + dy, f"{x:.2f}°", color=color)

# ------------ 4) Розрахунок ----------------------
FB = F_B(theta)
FC = F_C(theta)
cos_t = np.cos(theta)

H1 = FB
E1 = FB * cos_t

H2 = FB * FC
E2 = H2 * cos_t

zeros_H1 = find_zeros(H1)
zeros_E1 = zeros_H1[:]
zeros_H2 = find_zeros(H2)
zeros_E2 = find_zeros(E2)

max_H1 = find_maxima(H1)
max_E1 = find_maxima(E1)
max_H2 = find_maxima(H2)
max_E2 = find_maxima(E2)

lvl_H1 = level_cross_first_angle(H1)
lvl_E1 = level_cross_first_angle(E1)
lvl_H2 = level_cross_first_angle(H2)
lvl_E2 = level_cross_first_angle(E2)

# ============================================================
# ГРАФІК 1 — Однострижнева ДС (H и E)
# ============================================================
plt.figure(figsize=(10, 6))
ax = plt.gca()

H1n = norm(H1)
E1n = norm(E1)
cosn = norm(cos_t)

plt.plot(theta_deg, H1n, "k--", label="H-площина |F_B|")
plt.plot(theta_deg, cosn, "g-",  label="cosθ")
plt.plot(theta_deg, E1n, "b-",   label="E-площина |F_E|")
plt.axhline(HP_LEVEL, linestyle="--", color="gray", label="Рівень 0.707")

annotate_level(ax, lvl_H1, HP_LEVEL, color="black", dy=0.025)
annotate_level(ax, lvl_E1, HP_LEVEL, color="blue",  dy=-0.05)

annotate_zeros(ax, zeros_H1[:5], "θ0H", -0.10, "orange", "o", "black", n=5)
annotate_zeros(ax, zeros_E1[:5], "θ0E", -0.13, "purple", "s", "blue",  n=5)

annotate_maxima(ax, max_H1[:5], H1n, "θmH", "black", "x", n=5, dy=0.03)
annotate_maxima(ax, max_E1[:5], E1n, "θmE", "blue",  "*", n=5, dy=0.03)

plt.title("Однострижнева ДС (H та E)")
plt.xlabel("θ, град")
plt.ylabel("Нормована амплітуда")
plt.grid(True)
plt.ylim(0, 1.1)
plt.xlim(0, 90)
plt.legend()

# ============================================================
# ГРАФИК 2 — Двустрижнева ДС (H-площина)
# ============================================================
plt.figure(figsize=(10, 6))
ax = plt.gca()

H2n = norm(H2)
FBn = norm(FB)
FCn = norm(FC)

plt.plot(theta_deg, H2n, "b-",  label="|F_H(θ)| для 2 стрижнів")
plt.plot(theta_deg, FBn, "k--", label="|F_B(θ)| (один стрижень)")
plt.plot(theta_deg, FCn, "g-",  label="|F_C(θ)| (множник грат)")
plt.axhline(HP_LEVEL, linestyle="--", color="gray", label="Рівень 0.707")

annotate_level(ax, lvl_H2, HP_LEVEL, color="blue", dy=0.03)

annotate_zeros(ax, zeros_H2, "θ0H", -0.10, "orange", "o", "blue", n=10)
annotate_maxima(ax, max_H2, H2n, "θmH", "black", "x", n=10, dy=0.025)

plt.title("Двострижнева ДС — H-площина")
plt.xlabel("θ, град")
plt.ylabel("Нормована амплітуда")
plt.grid(True)
plt.ylim(0, 1.15)
plt.xlim(0, 90)
plt.legend()

# ============================================================
# ГРАФИК 3 — Двустрижнева ДС (E-площина)
# ============================================================
plt.figure(figsize=(10, 6))
ax = plt.gca()

E2n = norm(E2)
H2n = norm(H2)
FBn = norm(FB)
FCn = norm(FC)
cosn = norm(cos_t)

plt.plot(theta_deg, FBn, "k--", label="|F_B(θ)| (один стрижень)")
plt.plot(theta_deg, FCn, "g-",  label="|F_C(θ)| (множник грат)")
plt.plot(theta_deg, H2n, "b-",  label="|F_H(θ)| = |F_B·F_C|")
plt.plot(theta_deg, E2n, "r-", linewidth=2, label="|F_E(θ)| = |F_B·F_C·cosθ|")
plt.plot(theta_deg, cosn, "m-", label="|cos θ|")
plt.axhline(HP_LEVEL, linestyle="--", color="gray", label="Рівень 0.707")

annotate_level(ax, lvl_E2, HP_LEVEL, color="red", dy=0.03)

annotate_zeros(ax, zeros_E2, "θ0E", -0.07, "purple", "o", "red", n=10)
annotate_maxima(ax, max_E2, E2n, "θmE", "red", "*", n=10, dy=0.01)

# дублируем H2 нули/макс (как у тебя было)
annotate_zeros(ax, zeros_H2, "θ0H", -0.10, "orange", "o", "blue", n=10)
annotate_maxima(ax, max_H2, H2n, "θmH", "blue", "x", n=10, dy=0.03)

plt.title("Двострижнева ДС — E-площина")
plt.xlabel("θ, град")
plt.ylabel("Нормована амплітуда")
plt.grid(True)
plt.ylim(0, 1.15)
plt.xlim(0, 90)
plt.legend()

# ---------------------- Таблиці ----------------------
print_table("ТАБЛ. 1 — Нульові кути H₁", zeros_H1)
print_table("ТАБЛ. 2 — Нульові кути E₁", zeros_E1)
print_table("ТАБЛ. 3 — Максимуми H₁", max_H1, values=norm(H1))
print_table("ТАБЛ. 4 — Максимуми E₁", max_E1, values=norm(E1))

print_table("ТАБЛ. 5 — Нульові кути H₂", zeros_H2)
print_table("ТАБЛ. 6 — Нульові кути E₂", zeros_E2)
print_table("ТАБЛ. 7 — Максимуми H₂", max_H2, values=norm(H2))
print_table("ТАБЛ. 8 — Максимуми E₂", max_E2, values=norm(E2))

print("\nТАБЛ. 9 — Кути на рівні 0.707")
print("--------------------------------------")
print("| Графік           | θ (град)        |")
print("--------------------------------------")
if lvl_H1 is not None: print(f"| H₁                |  {lvl_H1:8.4f} |")
if lvl_E1 is not None: print(f"| E₁                |  {lvl_E1:8.4f} |")
if lvl_H2 is not None: print(f"| H₂                |  {lvl_H2:8.4f} |")
if lvl_E2 is not None: print(f"| E₂                |  {lvl_E2:8.4f} |")
print("--------------------------------------")

plt.show()