import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(0,2, 1001)

r2 = np.linspace(0, 2, 41)
theta2 = np.linspace(0,np.pi, 41)
Theta2, R2 = np.meshgrid(theta2, r2)
func_values = np.zeros_like(Theta2)

def f(r2,theta2, r1):
    nominator = np.exp(-4*(r1+r2))*r1**2*r2**2
    denominator = np.sqrt(r1**2 + r2**2 - 2*r1*r2*np.cos(theta2))
    if denominator == 0:
        return 0
    return nominator/denominator

F = np.vectorize(f)

for i in range(len(func_values.flat)):
    func_values.flat[i] = f(R2.flat[i], Theta2.flat[i],r1 = 1)

max_value = max(func_values.flat)
index_max = np.where(func_values == max_value)
optimal_choice_theta2 = Theta2[index_max]
optimal_choice_r2 = R2[index_max]
print("r2 = ", optimal_choice_r2)
print("theta2 = ", optimal_choice_theta2)

r1 = np.linspace(0.6,1.4,1001)

plt.plot(r1, F(1, optimal_choice_theta2, r1))
plt.xlabel("r")
plt.ylabel("wavefunction")
plt.show()
