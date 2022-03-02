# Import required libraries

import numpy as np

# Import amino acid sequence data

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

# Generation of A matrix

A = np.zeros((20,20))

for i in range(0,N):
  for j in range(i+1,N):
    for k in range(0,L):
      AA1 = AAdata[i][k]
      AA2 = AAdata[j][k]
      if AA1 != AA2:
        ind1 = AAtoIND[ord(AA1)]
        ind2 = AAtoIND[ord(AA2)]
        A[ind1][ind2] += 1
        A[ind2][ind1] += 1

np.set_printoptions(linewidth=100)
print(A)

# Calculation of the  λ  parameter with  λ  = 0.01 * Ntot/Atot

Atot = np.sum(np.sum(A)).astype(float)

Ntot = np.float(N * L)

lambd = 0.01 * Ntot / Atot

# Calculation of substitution matrix  M

Ncount = np.zeros(20,)

for i in range(0,20):
  for j in range(0,N):
    Ncount[i] += AAdata[j].count(AAletters[i])

print(Ncount) 
                               
M = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    if i !=j:
      M[i][j] = lambd * A[i][j] / Ncount[i]

M = M + np.diag(1 - np.sum(M,axis = 1))

np.set_printoptions(linewidth=1000)

#Sum of the rows of M matrix must be equal to 1. To check:
print(np.sum(M, axis = 1))

#print(M)
print(np.round(10000*M).astype(int))

# Higher order  M  matrices

M2 = np.matmul(M,M)
M3 = np.matmul(M2,M)
M4 = np.matmul(M2,M2)
M16 = np.matmul(M4,M4)
print(np.round(10000*M3).astype(int))

# Calculation of relative frequency matrix  R  and PAM3 log-odds matrix  S

#pi is equal to appearance of an amino acid / number of total amino acids
pi = np.zeros(20,)

for i in range(0,20):
  pi[i] = Ncount[i]/Ntot

print(pi)

# Relative frequency = (M^n) i,j/pi(j) 

R = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    R[i][j] = (M3[i][j] / pi[j])

print(np.round(100*R).astype(int))  


# Log-odds matrix S is  𝑆𝑖,𝑗 = 10 log10 𝑅𝑖,𝑗

S = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    S[i][j] = 10 * np.log10(R[i][j])

#print(S)
print("\n")
print(np.round(S).astype(int))
