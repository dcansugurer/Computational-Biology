# Import required libraries

import numpy as np

# Import input amino acid sequence data

filename = 'AAseq.txt'
file = open(filename,'r')
AAdata = file.read().splitlines()
file.close()
N = len(AAdata)
L = len(AAdata[0])
print("The file %s has %d sequences, each with %d letters for different amino acids" %(filename, N, L) )
print(L)
print(AAdata)

# Amino acid index matrix

#Calculate the A matrix

AAtoIND = np.zeros((100,),dtype=int)

AAletters = "ARNDCQEGHILKMFPSTWYV"

for i in range(0,20):
  AAtoIND[ord(AAletters[i])] = i

np.set_printoptions(linewidth=None)
print(AAtoIND)

# Calculation of  q  matrix by  q  =  Ai,j/Atot

# q i,j = A i,j/Atot
#A matrix for BLOSUM:

A_blos = np.zeros((20,20))

for i in range(0,N):
  for j in range(i+1,N):
    for k in range(0,L):
      AA1 = AAdata[i][k]
      AA2 = AAdata[j][k]
      if AA1 != AA2:
        ind1 = AAtoIND[ord(AA1)]
        ind2 = AAtoIND[ord(AA2)]
        A_blos[ind1][ind2] += 1
        A_blos[ind2][ind1] += 1
      else:
        A_blos[ind1][ind1] += 2

#print(A_blos)
print(np.round(A_blos).astype(int))

Atot_blos = np.sum(np.sum(A_blos)).astype(float)

q = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    q[i][j] = A_blos[i][j] / Atot_blos

print("\n")
print(np.round(q, 2).astype(float))

# Calculation of relative frequency matrix  R  by  Ri,j  =  qi,j/πiπj

R_blos = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    R_blos[i][j] = q[i][j] / (pi[i] * pi[j])

print(np.round(R_blos).astype(int))

#I converted zero values in R matrix to 1 to deal with -inf in S matrix
R2 = np.where(R_blos==0,1,R_blos)

#print(R_blos)  
#print(R2) 

# Calculation of log-odds matrix  S  for BLOSUM using  Si,j  = [2  log2   Ri,j ]

S_blos = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    S_blos[i][j] = 10 * np.log10(R2[i][j])

#print(S_blos)   
print(np.round(S_blos).astype(int))

