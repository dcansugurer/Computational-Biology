# import required libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sympy import symbols, solve
%matplotlib inline
fig = plt.figure(figsize=(10,4))

# parameters for positive autoregulation model
# k < beta / alpha
beta = 12
k = 14
n = 5
alpha = 0.4

x = np.linspace(0, 50)
prod_rate = beta * x**n / (x**n + k**n) # equeation for production rate of X
degr_rate = alpha * x # equation for degradation rate of X

# ploting of degradation and production rates of X
plt.plot(x, degr_rate, label='Degredation rate', color='r', zorder=2)
plt.plot(x, prod_rate, label='Production rate', color='g', zorder=2)
plt.title('Positive Autoregulation Model', fontsize=15)
plt.xlabel('[X] \n \n Graph 1: Kinetic model of positive autoregulation of X', fontsize=12)
plt.ylabel('Rate of Production & Degradation', fontsize=12)
plt.grid()
plt.xlim(-1, 40)
plt.ylim(-1, 15)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# calculation of critical [X] value
# solving the equation where production rate is equal to degradation rate
x = symbols('x')
equation = (alpha * (x ** (n+1))) - (beta * (x**n)) + (alpha * x * (k**n))
crit_value = solve(equation)
my_x = crit_value[2]
y = alpha * my_x
print(crit_value)
print("\n")
print("Critical [X] value for this positive autoregulation model is", round(my_x, 3))
print("\n")
plt.scatter(my_x, y, s=25, color='black', zorder=3)
plt.show()
