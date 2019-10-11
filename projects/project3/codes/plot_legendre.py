import matplotlib.pyplot as plt
import numpy as np


def func(x):
    return np.exp(-4*x)/x

r = np.linspace(0.05, 2, 100000)


plt.plot(r, func(r))
plt.show()
