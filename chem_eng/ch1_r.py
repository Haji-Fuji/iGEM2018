import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

np.set_printoptions(threshold=np.inf)

def function(V, c_A0, t, c_A, F):
	return V*diff(c_A0, t) = F*c_A0 - F*c_A

V = 10
F = 1
c_A0 = 1

t = np.arange(0, 50, 0.1, dtype=float) 



c_A = obeint(function,	t, args=(V, c_A0, F))

plt.plot(t, c_A)
plt.show()