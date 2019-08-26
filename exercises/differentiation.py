import sys
import numpy as np
import matplotlib.pyplot as plt

filename = str(sys.argv[1])
with open(filename, "r") as infile:
    first_line = infile.readline()
    h = []
    derivative = []
    error = []
    lines = infile.readlines()
    for line in lines:
        numbers = line.split()
        h.append(float(numbers[0]))
        derivative.append(float(numbers[1]))
        error.append(float(numbers[-1]))


h = np.array(h)
derivative = np.array(derivative)
error = np.array(error)

index = np.where(error == min(error[:20]))
print(h[index])


plt.plot(np.log10(h), error, label="Error as function of stepsize h")
plt.xlabel("log10(h)")
plt.ylabel("log10((f'(x) - f_approx'(x))/f'(x))")
plt.legend()
plt.show()
