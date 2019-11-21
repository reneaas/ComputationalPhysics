import matplotlib.pyplot as plt
import numpy as np


x = []
u = []

infilename = "test.txt"

with open(infilename, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        values = line.split()
        x.append(float(values[0]))
        u.append(float(values[1]))


def exact(x, t, N = 100):
    sum = 0
    for i in range(1,N+1):
        sum += ((-1)**i)/i * np.sin(i*np.pi*x)*np.exp(-(i*np.pi)**2 * t)
    sum *= 2/np.pi
    return sum + x

dx = 0.01

xx = np.linspace(dx,1-dx,1001)
t = 0.001

plt.plot(x,u, label = "approximation")
plt.plot(xx, exact(xx, t), "--", label = "exact")
plt.legend()
plt.show()
