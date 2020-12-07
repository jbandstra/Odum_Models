# import statements
import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# set global constants
E = 1
K1, K2, K3 = 0.07, 0.08, 0.002, 
K4, K5, K6 = 0.001, 0.002, 0.001
r1, r2 = K1*E, K2*E
Cap1, Cap2 = r1/K3, r2/K4
# set up the derivative function
def INTERACT(t, Q):
    Q1, Q2 = Q
    dQ1 = r1*Q1*(1 - Q1/Cap1) - K5*Q1*Q2
    dQ2 = r2*Q2*(1 - Q2/Cap2) - K6*Q1*Q2
    return [dQ1, dQ2]
# run the ODE solver
soln = solve_ivp(INTERACT, (0, 319), [3,3], dense_output=True)
# graph the populations vs time
t_arr = np.linspace(0,319,200)
plt.plot(t_arr, soln.sol(t_arr)[0], '-', label='Q_1')
plt.plot(t_arr, soln.sol(t_arr)[1], '--', label='Q_2')
plt.xlabel("Time")
plt.ylabel("Population Biomass")
plt.legend()
plt.figure()
plt.plot(soln.sol(t_arr)[0], soln.sol(t_arr)[1])
plt.xlabel("Q1 Biomass")
plt.ylabel("Q2 Biomass")
