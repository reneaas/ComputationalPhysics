import os
import sys
import matplotlib.pyplot as plt
import numpy as np

dataset = str(sys.argv[1])


#Empty lists to store computed data
MC_cycles = []                              #
E = []
E_sq = []
Mabs = []
Mabs_sq = []
M = []
M_sq = []
Cv = []
chi = []


if dataset == "L=2":
