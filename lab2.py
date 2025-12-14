import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


lambd = 3.2 / 100
a = 13 / 100

F1h = [1]
FC = [1]
FH = [1]
steps = [0]
SGP1 = 0
fS1 = 0
min_x_FH = []
min_y_FH = []
max_x_FH = []
max_y_FH = []

for theta in np.arange(0.01, np.pi / 2, 0.00001):
    mn1 = abs((1 + np.cos(theta)) / 2)
    mn2 = abs(math.cos((np.pi * a * math.sin(theta)) / lambd) / (1 - ((2 * a * math.sin(theta)) / lambd) ** 2))
    mn3 = mn1 * mn2
    F1h.append(mn1)
    FC.append(mn2)
    FH.append(mn3)

    if 0.707 < mn3 < 0.708:
        SGP1 = 2 * math.degrees(theta)
        fS1 = mn3
    steps.append(math.degrees(theta))

for i in range(1, len(FH) - 1):
    if FH[i] > FH[i - 1] and FH[i] > FH[i + 1]:
        max_x_FH.append(steps[i])
        max_y_FH.append(FH[i])
    if FH[i] < FH[i - 1] and FH[i] < FH[i + 1]:
        min_x_FH.append(steps[i])
        min_y_FH.append(0)

print("Табл. 1 - Значення нульових кутів")
print("------------------------------")
print("| № |  θ  |FH(θ)|")
for i, theta in enumerate(min_x_FH, start=1):
    print(f"| {i} | {theta:.2f} | 0 |")
print("------------------------------")

print("Табл. 2 - Значення максимальних кутів")
print("------------------------------")
print("| № |  θ  |FH(θ)|")
for i, theta in enumerate(max_x_FH, start=1):
    print(f"| {i} | {theta:.2f} | {max_y_FH[i-1]:.2f} |")
print("------------------------------")

print("Ширина головної пелюстки в площині H = " + str(round(SGP1, 2)) + '\u00b0')

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(steps, F1h, linewidth=0.7, label="$ F_{1h}(\theta) $")
ax.plot(steps, FC, linewidth=0.7, label="$ F_{C}(\theta) $")
ax.plot(steps, FH, linewidth=0.7, label="$ F_{H}(\theta) $")
ax.plot(SGP1/2, fS1, 'ro', markersize=4, label="ШГП в площині H")
ax.annotate(f'({fS1:.3f}, {SGP1/2:.2f}\u00b0)',
            xy=(SGP1 / 2, fS1),
            xytext=((SGP1/2)+3, fS1+0.05),
            arrowprops=dict(arrowstyle='->', color='black'), fontsize=8)
ax.plot(max_x_FH, max_y_FH, "o", markersize=4, color="black", label="$ \theta_{max} $ H")
ax.plot(min_x_FH, min_y_FH, "o", markersize=4, color="blue", label="$ \theta_{min} $ H")
ax.set_xlabel('θ' + '\u00b0', fontsize=10)
ax.set_ylabel('|Fh(θ' + '\u00b0' + ')|, |F1h(θ' + '\u00b0' + ')|, |Fc(θ' + '\u00b0' + ')|', fontsize=10)
plt.xticks(np.arange(0, 100, 2), fontsize=7)
plt.yticks(np.arange(0, 1.2, 0.1), fontsize=7)
plt.ylim(-0.01, 1.01)
plt.xlim(0, 90.5)
plt.legend(loc="upper right", fontsize=7)
plt.grid(which='both', linestyle='--', linewidth=0.2, color='gray')
plt.title(f'ДС пірамідального рупора в площині Н: λ = {lambd*100} см, a = {a*100} см', fontsize=10)
fig.savefig("ДС_рупора_варіант_8.jpg", dpi=600)
plt.show()