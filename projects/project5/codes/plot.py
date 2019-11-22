import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def exact(x, t, N = 1000):
    sum = 0
    for i in range(1,N+1):
        sum += ((-1)**i)/i * np.sin(i*np.pi*x)*np.exp(-(i*np.pi)**2 * t)
    sum *= 2/np.pi
    return sum + x


dx = 0.01
xx = np.linspace(dx,1-dx,1001)
t = 0.005

methods = ["explicit", "implicit", "CN"]

for method in methods:
    x = []
    u = []
    infilename = "test_" + str(method) + ".txt"
    with open(infilename, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            x.append(float(values[0]))
            u.append(float(values[1]))

    plt.plot(x,u, label = str(method))
plt.plot(xx, exact(xx, t), "--", label = "exact")
plt.legend()
plt.show()
