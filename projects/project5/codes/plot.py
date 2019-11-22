import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from mpl_toolkits import mplot3d

d = int(input("d?"))

if d == 1:
    def exact(x, t, N = 1000):
        sum = 0
        for i in range(1,N+1):
            sum += ((-1)**i)/i * np.sin(i*np.pi*x)*np.exp(-(i*np.pi)**2 * t)
        sum *= 2/np.pi
        return sum + x

    dx = float(input("Give dx: [0.1 or 0.01] \n"))
    time = int(input("Initial time [type 1] or steady state [type 2]: \n"))
    xx = np.linspace(dx,1-dx,1001)


    methods = ["explicit", "implicit", "CN"]
    path = "results/1D/"


    for method in methods:
        x = []
        u = []
        infilename = str(method) + "_dx_" + str(dx) + "_time_" + str(time) + ".txt"
        with open(path + infilename, "r") as infile:
            t = float(infile.readline().split()[0])
            lines = infile.readlines()
            for line in lines:
                values = line.split()
                x.append(float(values[0]))
                u.append(float(values[1]))

        plt.plot(x,u, label = str(method))
    plt.plot(xx, exact(xx, t), "--", label = "exact")
    plt.legend()
    plt.show()

if d == 2:
    x = [(i+1)*0.1 for i in range(98)]
    y = [(i+1)*0.1 for i in range(98)]
    u = np.zeros((len(x), len(y)))
    infilename = "results/2D/2D_Results.txt"
    with open(infilename, "r") as infile:
        t = float(infile.readline().split()[0])
        lines = infile.readlines()
        for i in range(len(lines)):
            values = lines[i].split()
            for j in range(len(values)):
                u[i,j] = float(values[j])
    u = u.T;


    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf = ax.contour3D(x, y, u, 98)
    ax.set_zlim(0, 1)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
