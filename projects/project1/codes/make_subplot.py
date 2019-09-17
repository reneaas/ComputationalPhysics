import numpy as np
import matplotlib.pyplot as plt
import sys
from os import system

image_address = ["results/general_algorithm/plots/max_error_general.png", \
                "results/special_algorithm/plots/max_error_special.png",\
                "results/LU/plots/max_error_LU.png"]

images = [plt.imread(i) for i in image_address]

fig = plt.figure()

for i in range(len(images)-1):
    fig.add_subplot(1,2,i+1)
    plt.imshow(images[i])
    plt.axis("off")


figurename = "subplot_error.png"
plt.savefig(figurename, dpi = 1000)
plt.close()

directory = "~/Documents/skole/comphys/projects/project1/codes/results"
system("mv" + " " + figurename + " " + directory)
