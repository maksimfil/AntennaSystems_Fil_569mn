import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
F = 560e6  # Частота в Гц
c = 3e8  # Швидкість світла в м/с
N = 9  # Кількість елементів

lambda_ = c / F
print(f"# Розрахунок робочої довжини хвилі = {lambda_}")

Dp_min, Dp_max = 0.1 * lambda_, 0.25 * lambda_
Dd_min, Dd_max = 0.1 * lambda_, 0.34 * lambda_

theta = np.linspace(-90, 90, 360)  # Кути в градусах для кращого аналізу

theta_rad = np.radians(theta)

F1H_theta = np.ones_like(theta_rad)
Fc_theta = np.sinc(N * np.sin(theta_rad))
FH_theta = F1H_theta * Fc_theta

F1E_theta = np.cos((np.pi / 2) * np.sin(theta_rad)) / np.cos(theta_rad)
FE_theta = F1E_theta * FH_theta

FH_theta /= np.max(FH_theta)
FE_theta /= np.max(FE_theta)

threshold = 0.707

def beamwidth(theta, F_theta):
    indices = np.where(F_theta >= threshold)[0]
    if len(indices) > 1:
        return theta[indices[-1]] - theta[indices[0]]
    return None

beamwidth_H = beamwidth(theta, FH_theta)
beamwidth_E = beamwidth(theta, FE_theta)

sidelobe_H = np.max(FH_theta[np.where(FH_theta < 1)])
sidelobe_E = np.max(FE_theta[np.where(FE_theta < 1)])

plt.figure(figsize=(8, 6))
plt.polar(np.radians(theta), FH_theta, label=f'H-площина (Ширина ГП: {beamwidth_H:.2f}°)')
plt.polar(np.radians(theta), FE_theta, linestyle='dashed', label=f'E-площина (Ширина ГП: {beamwidth_E:.2f}°)')
plt.legend()
plt.title(f'Діаграма спрямованості\nРівень бічної пелюстки H: {sidelobe_H:.2f}, E: {sidelobe_E:.2f}')
plt.show()
