# Import required libraries

import numpy as np
import scipy as sp
from scipy import linalg
import matplotlib.pyplot as plt
import math
from math import e

# Random transition matrix generation

rand_rates = np.random.rand(4,4) * (1 - np.eye(4))
rand_Q = rand_rates - np.diag(np.sum(rand_rates,1))
print(rand_Q)

# Calculate the transition probability matrix 𝑃(𝑡)corresponding to the transition rate matrix 𝑄 above 
# Plot the corresponding transition probabilities between each nucleotide pair in a 4×4table for 𝑡∈[0,5].

P = np.zeros((4,4,101))

ind=0;
for t in np.arange(0,10.1,0.1):
  cP = sp.linalg.expm(rand_Q*t)
  P[:,:,ind] = cP
  ind+=1

print(cP)

# Plotting of probability curves as functions of time

fig = plt.figure(figsize=(20,10))

ind = 0
for (i,j) in [(i,j) for i in range(0,4) for j in range(0,4)]:
  ind +=1 
  axs = fig.add_subplot(4,4,ind)
  axs.plot(np.linspace(0,10,101), P[i,j,:])
  axs.legend(['P(%d,%d)' % (i+1, j+1)])
  axs.set_xlim([0,5])
  axs.set_ylim([0,1])
  axs.grid()
  if i == 3:
    axs.set_xlabel('t')
  else: 
      axs.set_xticklabels([])
      
# Extention of time span to 𝑡∈[0,100] to see convergence of the probabilty curves

fig_2 = plt.figure(figsize=(20,10))

ind = 0
for (i,j) in [(i,j) for i in range(0,4) for j in range(0,4)]:
  ind +=1 
  axs = fig_2.add_subplot(4,4,ind)
  axs.plot(np.linspace(0,100,101), P[i,j,:])
  axs.legend(['P(%d,%d)' % (i+1, j+1)])
  axs.set_xlim([0,100])
  axs.set_ylim([0,1])
  axs.grid()
  if i == 3:
    axs.set_xlabel('t')
  else: 
      axs.set_xticklabels([])
      
# Functional relationship between the sequence distance 𝐷 and the evolutionary distance 𝑑 between two nucleotide sequences.
# Relationship as a graph of 𝑑 versus 𝐷.

# at refers to alpha*t 
A = np.zeros((4,4,101))
ind=0;
for t in np.arange(0,10.1,0.1):
  at = rand_Q*t
  A[:,:,ind] = at
  ind+=1

# total evolutionary distance: d(t) = 6𝛼𝑡
d = 6*A

#Sequence distance through the 2t evolutionary time : D(t) = 3/4 - 3/4*e**(-8𝛼𝑡)
D = -0.75*e**(-8*A) + 0.75


# Axis limits of the plot are set for y:[0,5], x:[0,1] to have similar plot with slide 18
# Diagonal plots have negative values, we cannot see plots in that setting.
fig_D = plt.figure(figsize=(20,10))

ind = 0
for (i,j) in [(i,j) for i in range(0,4) for j in range(0,4)]:
  ind +=1 
  axs = fig_D.add_subplot(4,4,ind)
  axs.plot(D[i,j], d[i,j])
  axs.legend(['D(%d,%d)' % (i+1, j+1)])
  axs.set_ylim([0,5])
  axs.set_xlim([0,1])
  axs.grid()
  axs.set_ylabel('d')
  if i == 3:
    axs.set_xlabel('D')
  else: 
    axs.set_xticklabels([]) 
    
# at refers to alpha*t 
A = np.zeros((4,4,101))
ind=0;
for t in np.arange(0,10.1,0.1):
  at = rand_Q*t
  A[:,:,ind] = at
  ind+=1

# total evolutionary distance: d(t) = 6𝛼𝑡
d = 6*A

# Sequence distance through the 2t evolutionary time : D(t) = 3/4 - 3/4*e**(-8𝛼𝑡)
D = -0.75*e**(-8*A) + 0.75

# Plots without axis limit
fig_D = plt.figure(figsize=(20,10))

ind = 0
for (i,j) in [(i,j) for i in range(0,4) for j in range(0,4)]:
  ind +=1 
  axs = fig_D.add_subplot(4,4,ind)
  axs.plot(D[i,j], d[i,j])
  axs.legend(['D(%d,%d)' % (i+1, j+1)])
  axs.grid()
  axs.set_ylabel('d')
  if i == 3:
    axs.set_xlabel('D')
  else: 
      axs.set_xticklabels([])
