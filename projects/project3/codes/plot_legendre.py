import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(1, 3, 101)
eps = 0.1

def f(r):
    return np.exp(-4*r)/r


F = f(r)
indices = np.where(f(r) < f(2.7))
func_value = F[indices]
print(func_value)


#plt.plot(r, f(r))


#plt.plot(r, f(r), label = "f(r)")
plt.plot(r[indices], func_value)
plt.show()
