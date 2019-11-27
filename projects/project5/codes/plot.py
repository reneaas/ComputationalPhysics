import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from mpl_toolkits import mplot3d
plt.rc("text", usetex = True)

d = int(input("d? \n"))

if d == 1:

    def exact_1D(x, t, N = 1000):
        """
        Analytical solution in the 1D-case.
        """
        sum = 0
        for i in range(1,N+1):
            sum += ((-1)**i)/i * np.sin(i*np.pi*x)*np.exp(-(i*np.pi)**2 * t)
        sum *= 2/np.pi
        return sum + x



    fig1 = plt.figure(); figurename1 = "dx_0.1_time_1.pdf"
    fig2 = plt.figure(); figurename2 = "dx_0.1_time_2.pdf"
    fig3 = plt.figure(); figurename3 = "dx_0.01_time_1.pdf"
    fig4 = plt.figure(); figurename4 = "dx_0.01_time_2.pdf"

    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax3 = fig3.add_subplot(111)
    ax4 = fig4.add_subplot(111)


    methods = ["explicit", "implicit", "CN"]
    path = "results/1D/"

    dx = [0.1, 0.01]

    dx1 = {}
    dx2 = {}
    t = []

    for i in range(len(dx)):
        for j in range(1,3):
            for method in methods:
                x = []
                u = []
                infilename = method + "_dx_" + str(dx[i]) + "_time_" + str(j) + ".txt"
                with open(path + infilename, "r") as infile:
                    t.append(float(infile.readline().split()[0]))
                    lines = infile.readlines()
                    for line in lines:
                        values = line.split()
                        x.append(float(values[0]))
                        u.append(float(values[1]))
                if i == 0:
                    dx1["x_" + method + "_time_{0}".format(j)] = x
                    dx1["u_" + method + "_time_{0}".format(j)] = u
                else:
                    dx2["x_" + method + "_time_{0}".format(j)] = x
                    dx2["u_" + method + "_time_{0}".format(j)] = u

    time = t[::3]
    x1 = np.linspace(dx[0],1-dx[0],1001)
    x2 = np.linspace(dx[1],1-dx[1],1001)



    for i in range(3):
        ax1.plot(dx1["x_" + methods[i] + "_time_1"], dx1["u_" + methods[i] + "_time_1"], label = methods[i])
        ax1.set_xlabel("$x$", size = 14)
        ax1.set_ylabel("$u(x,t)$", size = 14)
        ax1.tick_params(labelsize = 15)

        ax2.plot(dx1["x_" + methods[i] + "_time_2"], dx1["u_" + methods[i] + "_time_2"], label = methods[i])
        ax2.set_xlabel("$x$", size = 14)
        ax2.set_ylabel("$u(x,t)$", size = 14)
        ax2.tick_params(labelsize = 15)

        ax3.plot(dx2["x_" + methods[i] + "_time_1"], dx2["u_" + methods[i] + "_time_1"], label = methods[i])
        ax3.set_xlabel("$x$", size = 14)
        ax3.set_ylabel("$u(x,t)$", size = 14)
        ax3.tick_params(labelsize = 15)

        ax4.plot(dx2["x_" + methods[i] + "_time_2"], dx2["u_" + methods[i] + "_time_2"], label = methods[i])
        ax4.set_xlabel("$x$", size = 14)
        ax4.set_ylabel("$u(x,t)$", size = 14)
        ax4.tick_params(labelsize = 15)

    ax1.plot(x1, exact_1D(x1,time[0]), "--", label = "exact")
    ax1.legend(fontsize = 12)
    ax1.set_title("Time = {0}, r = 0.5".format(time[0]))
    fig1.savefig(figurename1)

    ax2.plot(x1, exact_1D(x1,time[1]), "--", label = "exact")
    ax2.legend(fontsize = 12)
    ax2.set_title("Time = {0}, r = 0.5".format(time[1]))
    fig2.savefig(figurename2)

    ax3.plot(x2, exact_1D(x2,time[2]), "--", label = "exact")
    ax3.legend(fontsize = 12)
    ax3.set_title("Time = {0}, r = 0.5".format(time[2]))
    fig3.savefig(figurename3)

    ax4.plot(x2, exact_1D(x2,time[3]), "--", label = "exact")
    ax4.legend(fontsize = 12)
    ax4.set_title("Time = {0}, r = 0.5".format(time[3]))
    fig4.savefig(figurename4)


    os.system("mv" + " " + figurename1 + " " + figurename2 + " " + figurename3 + " " + figurename4 + " " + path)


if d == 2:
    def exact_2D(x, y, t, terms = 1001, a = 0, b = 1):
        """
        Analytic solution in the 2D-case.
        """
        s = 0
        for n in range(1,terms+1):
            for m in range(1,terms+1):
                coeff = 4*(a-b)*((-1)**(m+n) - (-1)**m)/(m*n*np.pi**2)\
                        - 4*a/(n*m*np.pi**2)*(1 - (-1)**m)*(1 - (-1)**n);
                s += coeff*np.sin(n*np.pi*x)*np.sin(m*np.pi*y)*np.exp(-t*(n**2 + m**2)*np.pi**2)
        return s + (b-a)*y + a;


    infilename = "results/2D/2D_Results.txt"
    with open(infilename, "r") as infile:
        t = float(infile.readline().split()[0])
        lines = infile.readlines()
        gridpoints = len(lines)
        u = np.zeros((gridpoints, gridpoints))
        for i in range(len(lines)):
            values = lines[i].split()
            for j in range(len(values)):
                u[i,j] = float(values[j])


    u = u.T                     #Transposes the data to order the matrix correctly for plotting.
    x = [(i+1)*0.01 for i in range(gridpoints)]
    y = [(i+1)*0.01 for i in range(gridpoints)]


    figurename1 = "numerical_2D.pdf"
    figurename2 = "analytical_2D.pdf"

    """
    #Plots the numerical solution z = u(x,y,t) as a surface plot in 3D for a specific time t.
    X,Y = np.meshgrid(x,y)
    ax1 = plt.axes(projection='3d')
    surf = ax1.contour3D(x, y, u, 98)
    ax1.set_zlim(0, 1)
    ax1.set_title("Time = {0}, r = 0.25".format(t))
    ax1.set_xlabel("x", size = 14)
    ax1.set_ylabel("y", size = 14)
    ax1.set_zlabel(r"$u(x,y,t_0)$", size = 14)
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    fig1.savefig(figurename1)
    plt.close()




    #Plots the analytical solution z = u(x,y,t) as a surface plot in 3D for a specific time t.
    xx = np.linspace(0,1,41)
    yy = np.linspace(0,1,41)
    X,Y = np.meshgrid(xx,yy)
    Z = exact_2D(X, Y,t = 0.0075, terms = 100)
    ax2 = plt.axes(projection='3d')
    surf = ax2.contour3D(X, Y, Z, 98)
    ax2.set_zlim(0, 1)
    ax2.set_title("Time = {0}, r = 0.25".format(t))
    ax2.set_xlabel("x", size = 14)
    ax2.set_ylabel("y", size = 14)
    ax2.set_zlabel(r"$u(x,y,t_0)$", size = 14)
    fig2.colorbar(surf, shrink=0.5, aspect=5)
    fig2.savefig(figurename2)
    plt.close()
    """

    #Plots the numerical solution z = u(x,y,t) as a surface plot in 3D for a specific time t.
    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf = ax.contour3D(x, y, u, 98)
    ax.set_zlim(0, 1)
    ax.set_title("Time = {0}, r = 0.25".format(t))
    ax.set_xlabel("$x$", size = 14)
    ax.set_ylabel("$y$", size = 14)
    ax.set_zlabel(r"$u(x,y,t_0)$", size = 14)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig(figurename1)
    plt.close()

    #Plots the analytical solution z = u(x,y,t) as a surface plot in 3D for a specific time t.
    xx = np.linspace(0,1,41)
    yy = np.linspace(0,1,41)
    X,Y = np.meshgrid(xx,yy)
    Z = exact_2D(X, Y,t = 0.0075, terms = 100)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf = ax.contour3D(X, Y, Z, 98)
    ax.set_zlim(0, 1)
    ax.set_title("Time = {0}, r = 0.25".format(t))
    ax.set_xlabel("$x$", size = 14)
    ax.set_ylabel("$y$", size = 14)
    ax.set_zlabel(r"$u(x,y,t_0)$", size = 14)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig(figurename2)
    plt.close()


    #Plots a contour plot of the numerical solution z = u(x,y,t) for a specific time t.
    figurename3 = "contour_2D.pdf"

    plt.contourf(x,y,u)
    plt.xlabel("$x$", fontsize = 14)
    plt.title(("Time = {0}, r = 0.25".format(t)), fontsize = 14)
    plt.ylabel("$y$", fontsize = 14)
    plt.colorbar().set_label("$u(x,y,t_0)$", size = 14)
    plt.savefig(figurename3)


    path = "results/2D/"
    os.system("mv" + " " + figurename1 + " " + figurename2 + " " + figurename3 + " " + path)
