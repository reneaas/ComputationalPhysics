import numpy as np



def avg_speedup(t1,t2):
    #t1, t2 are vectors
    N = len(t1)
    time = np.zeros(N)
    for i in range(N):
            time[i] = t1[i]/t2[i]

    return np.sum(time)/N


t1 = np.array([1,2,3])
t2 = np.array([4,5,5])

print(avg_speedup(t1,t2))
