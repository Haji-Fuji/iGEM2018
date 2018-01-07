from sympy import*
import matplotlib.pyplot as plt
import numpy as np

V = 10
F = 1
cA0 = 1

dsolve(V*diff(cA, t) = F*cA0 - F*cA)

t = np.arange(0, 50, 100)

plt.plot(t, cA)