import numpy as np
import matplotlib.pyplot as plt

infilename = "buckling_beam.txt"

x = []
u = []

def f(x, n=1):
    return np.sin(np.pi*n*x)/10



with open(infilename, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        vals = line.split()
        x.append(float(vals[0]))
        u.append(float(vals[1]))

N = len(x)
analytical = np.zeros(N)
for i in range(N):
    analytical[i] = np.sin(i*np.pi/N)

analytical = analytical/np.linalg.norm(analytical)
x = np.array(x)
#u = np.array(u)/np.linalg.norm(u)
plt.plot(x,u)
plt.plot(x,analytical, label= "analytical")
plt.legend()
plt.show()

ratio = analytical/u

print(np.mean(ratio))
