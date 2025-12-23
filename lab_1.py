import math
import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# ВХІДНІ ДАНІ
# Варіант №9: F = 580 МГц, N = 8
# ==========================================================
N = 8
F = 580e6
C = 299_792_458

# Параметри обчислення
THETA_START = 0.01
THETA_STOP = np.pi / 2
THETA_STEP = 1e-5
HP_MIN = 0.707
HP_MAX = 0.708


def main():
    # Хвильові параметри
    lambd = C / F
    d = 0.25 * lambd
    k = 2 * np.pi / lambd

    print(f"λ (довжина хвилі) = {lambd:.4f} (m)")
    print(f"dcp (крок елементів) = {d:.4f} (m)")
    print(f"k (хвильове число) = {k:.4f} (rad/m)")

    # Дані для графіків
    F1E = [1.0]
    FC = [1.0]
    FE = [1.0]
    steps = [0.0]

    # ШГП на рівні 0.707
    SGP1 = 0.0
    SGP2 = 0.0
    fS1 = 0.0
    fS2 = 0.0

    # Екстремуми
    max_x_FC, max_y_FC = [], []
    max_x_FE, max_y_FE = [], []
    min_x_FC, min_y_FC = [], []
    min_x_FE, min_y_FE = [], []

    # ==========================================================
    # РОЗРАХУНОК ДС
    # ==========================================================
    for theta in np.arange(THETA_START, THETA_STOP, THETA_STEP):
        # F1e(theta) - множник елемента (площина E, диполь)
        mn1 = abs(np.cos((np.pi / 2) * np.sin(theta)) / np.cos(theta))

        # FH(theta) - множник решітки (площина H)
        arg_H = (k * d * (1 - np.cos(theta))) / 2
        mn2 = abs(np.sin(N * arg_H) / (N * np.sin(arg_H)))

        # FE(theta) - загальна ДС (площина E)
        mn3 = mn1 * mn2

        F1E.append(mn1)
        FC.append(mn2)
        FE.append(mn3)

        theta_deg = math.degrees(theta)
        steps.append(theta_deg)

        # Перше попадання в діапазон 0.707..0.708
        if HP_MIN < mn2 < HP_MAX and SGP1 == 0:
            SGP1 = 2 * theta_deg
            fS1 = mn2

        if HP_MIN < mn3 < HP_MAX and SGP2 == 0:
            SGP2 = 2 * theta_deg
            fS2 = mn3

    # ==========================================================
    # ЕКСТРЕМУМИ (бічні пелюстки та нулі)
    # ==========================================================
    for i in range(1, len(FC) - 1):
        # Максимуми
        if FC[i] > FC[i - 1] and FC[i] > FC[i + 1]:
            max_x_FC.append(steps[i])
            max_y_FC.append(FC[i])

        if FE[i] > FE[i - 1] and FE[i] > FE[i + 1]:
            max_x_FE.append(steps[i])
            max_y_FE.append(FE[i])

        # Мінімуми
        if FC[i] < FC[i - 1] and FC[i] < FC[i + 1]:
            min_x_FC.append(steps[i])
            min_y_FC.append(FC[i])

        if FE[i] < FE[i - 1] and FE[i] < FE[i + 1]:
            min_x_FE.append(steps[i])
            min_y_FE.append(FC[i])

    print("-----------------------------------")
    print("Ширина головної пелюстки (ШГП) на рівні 0.707:")
    print(f"  В площині H = {SGP1:.2f}°")
    print(f"  В площині E = {SGP2:.2f}°")
    print("-----------------------------------")

    # ==========================================================
    # ТАБЛИЦЯ
    # ==========================================================
    max_len = max(len(max_x_FC), len(max_x_FE), len(min_x_FC), len(min_x_FE))
    print("Табл. 1 - Аналіз ДС Директорної антени в площині Н та Е")
    print("---------------------------------------------------------")
    print("| № |θ_min_H|θ_min_E|θ_max_H| F_H(θ) |θ_max_E| F_E(θ) |")
    print("---------------------------------------------------------")
    for idx in range(max_len):
        v1 = f"{idx + 1}"
        v2 = f"{min_x_FC[idx]:.2f}" if idx < len(min_x_FC) else "  -  "
        v3 = f"{min_x_FE[idx]:.2f}" if idx < len(min_x_FE) else "  -  "
        v4 = f"{max_x_FC[idx]:.2f}" if idx < len(max_x_FC) else "  -  "
        h4 = f"{max_y_FC[idx]:.3f}" if idx < len(max_y_FC) else "  -  "
        v5 = f"{max_x_FE[idx]:.2f}" if idx < len(max_x_FE) else "  -  "
        h5 = f"{max_y_FE[idx]:.3f}" if idx < len(max_y_FE) else "  -  "
        print(f"|{v1:^2}|{v2:^7}|{v3:^7}|{v4:^7}|{h4:^8}|{v5:^7}|{h5:^8}|")
    print("---------------------------------------------------------")

    # ==========================================================
    # ГРАФІК
    # ==========================================================
    fig, ax = plt.subplots(figsize=(20 / 2.54, 12 / 2.54))

    ax.plot(steps, F1E, linewidth=0.7, label="$ F_{1e}(\\theta) $ - Множник елемента")
    ax.plot(steps, FC, linewidth=0.7, label="$ F_{H}(\\theta) $ - Площина H (Множник решітки)")
    ax.plot(steps, FE, linewidth=0.7, label="$ F_{E}(\\theta) $ - Площина E (Загальна ДС)")

    # Позначки ШГП
    ax.plot(SGP1 / 2, fS1, "ro", markersize=4, label=f"ШГП в H ({SGP1:.2f}°)")
    ax.plot(SGP2 / 2, fS2, "go", markersize=4, label=f"ШГП в E ({SGP2:.2f}°)")

    # Анотації
    ax.annotate(
        f"({SGP1 / 2:.2f}°)",
        xy=(SGP1 / 2, fS1),
        xytext=(SGP1 / 2 + 3, fS1 - 0.1),
        arrowprops=dict(arrowstyle="->", color="black"),
        fontsize=6,
    )
    ax.annotate(
        f"({SGP2 / 2:.2f}°)",
        xy=(SGP2 / 2, fS2),
        xytext=(SGP2 / 2 - 12, fS2 + 0.05),
        arrowprops=dict(arrowstyle="->", color="black"),
        fontsize=6,
    )

    # Пунктирні лінії
    ax.hlines(y=fS1, xmin=0, xmax=SGP1 / 2, colors="r", linestyles="--", linewidth=0.5)
    ax.vlines(x=SGP1 / 2, ymin=0, ymax=fS1, colors="r", linestyles="--", linewidth=0.5)
    ax.vlines(x=SGP2 / 2, ymin=0, ymax=fS2, colors="g", linestyles="--", linewidth=0.5)

    # Екстремуми
    ax.plot(max_x_FC, max_y_FC, "o", markersize=4, color="black", label="$\\theta_{max}$ FH")
    ax.plot(max_x_FE, max_y_FE, "o", markersize=4, color="grey", label="$\\theta_{max}$ FE")
    ax.plot(min_x_FC, min_y_FC, "o", markersize=4, color="blue", label="$\\theta_{min}$")

    ax.set_xlabel("Кут відносно осі антени, $\\theta$ (°)", fontsize=10)
    ax.set_ylabel("Нормована ДС, $|F(\\theta)|$", fontsize=10)
    ax.set_title(f"Нормовані ДС Директорної антени", fontsize=12)

    plt.xticks(np.arange(0, 100, 5), fontsize=7)
    plt.yticks(np.arange(0, 1.2, 0.1), fontsize=7)
    plt.ylim(-0.01, 1.05)
    plt.xlim(0, 90.5)

    plt.legend(loc="upper right", fontsize=7)
    plt.grid(which="both", linestyle="--", linewidth=0.2, color="gray")

    fig.savefig(f"ДС_Director_Var9_N{N}_F580.jpg", dpi=600)
    plt.show()


if __name__ == "__main__":
    main()