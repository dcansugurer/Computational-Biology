# Import required libraries

import pandas as pd
import numpy as np

# Removing the gaps from alignment to generate PAM and BLOSUM matrices

filename = 'aligned.txt' #contains my multiple sequence alignment result
file = open(filename,'r')
aln = file.read().splitlines()
file.close()

Num_aln = len(aln)
Len_aln = len(aln[0])

temp = [] #list for storing column indexes contains "-" for gap regions

for i in range(0, Num_aln):
  for j in range(0, Len_aln):
    AA = aln[i][j]
    if AA == "-":
      temp.append(j)

temp.sort() #sort column indexes from smallest to largest


op = [] #removes duplicates and have only unique column indexes
for x in temp:
    if x not in op:
        op.append(x)

df = pd.DataFrame(aln)
my_df = df[0].apply(lambda x: pd.Series(list(x))) 

no_gap_df = my_df.drop((op), axis=1) #remove the columns that contain "-"

no_gap_df.to_csv('data_file.csv')
#After this point, NaN values were removed manually from data_file.csv and data with tab delimiter was converted into non-delimiter text file by using online tools.

filename = 'aligned_no_gap.txt'
file = open(filename,'r')
AAdata = file.read().splitlines()
file.close()
N = len(AAdata)
L = len(AAdata[0])
print("The file %s has %d sequences, each with %d letters for different amino acids" %(filename, N, L) )

#Calculate the A matrix
AAtoIND = np.zeros((100,),dtype=int)
AAletters = "ARNDCQEGHILKMFPSTWYV"

for i in range(0,20):
  AAtoIND[ord(AAletters[i])] = i

np.set_printoptions(linewidth=None)

# PAM
# Generate A matrix

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

np.set_printoptions(linewidth=1000)
#print(np.round(A).astype(int))

#  Calculation of the  λ  parameter with  λ  = 0.01  Ntot/Atot

Atot = np.sum(np.sum(A)).astype(float)

Ntot = np.float(N * L)

lambd = 0.01 * Ntot / Atot

# Calculation of substitution matrix  M

Ncount = np.zeros(20,)

for i in range(0,20):
  for j in range(0,N):
    Ncount[i] += AAdata[j].count(AAletters[i])
                               
M = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    if i !=j:
      M[i][j] = lambd * A[i][j] / Ncount[i]

M = M + np.diag(1 - np.sum(M,axis = 1))

np.set_printoptions(linewidth=1000)

#Sum of the rows of M matrix must be equal to 1. To check:
#print(np.sum(M, axis = 1))

#print(np.round(10000*M).astype(int))

# Higher order M matrices

#Calculation of higher order M matrices to get rid of zeros in relative frequency matrix R

M2 = np.matmul(M,M)
M3 = np.matmul(M2,M)
#print(np.round(10000*M3).astype(int))

# Calculation of relative frequency matrix  R  and PAM3 log-odds matrix  S

#pi is equal to appearance of an amino acid / number of total amino acids
pi = np.zeros(20,)

for i in range(0,20):
  pi[i] = Ncount[i]/Ntot

#print(pi)

# Relative frequency = (M^n) i,j/pi(j) 

#R = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    R[i][j] = (M3[i][j] / pi[j])

#print(np.round(100*R).astype(int))  


# Log-odds matrix S is  𝑆𝑖,𝑗 = 10 log10 𝑅𝑖,𝑗

S = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    S[i][j] = 10 * np.log10(R[i][j])

#print(S)
#print(np.round(S).astype(int))

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 200)
pd.set_option('expand_frame_repr', False)

df_S_PAM = pd.DataFrame(np.round(S).astype(int), columns = ["A-0","R-1","N-2","D-3","C-4","Q-5","E-6","G-7","H-8","I-9","L-10","K-11","M-12","F-13","P-14","S-15","T-16","W-17","Y-18","V-19"], 
                        index = ["A-0","R-1","N-2","D-3","C-4","Q-5","E-6","G-7","H-8","I-9","L-10","K-11","M-12","F-13","P-14","S-15","T-16","W-17","Y-18","V-19"])
print("\n")
print("Log-odds matrix for PAM:")
print("\n")
print(df_S_PAM)
print("\n")


#masking diagonal and non-diagonal elements
mask_diag = np.ones((20,20)) 
mask_diag = (np.diag(np.ones(20))).astype(np.bool) 
mask_non_diag = np.ones((20,20)) 
mask_non_diag = (mask_non_diag - np.diag(np.ones(20))).astype(np.bool)

most_cons_PAM = np.where(S == np.amax(S[mask_diag]))
print("Most conserved aa index is:", most_cons_PAM)
least_cons_PAM = np.where(S == np.amin(S[mask_diag]))
print("Least conserved aa index is:", least_cons_PAM)

most_likely_subs = np.where(S == np.amax(S[mask_non_diag]))
print("Most likely substitution index is:", most_likely_subs)
least_likely_subs = np.where(S == np.amin(S[mask_non_diag]))
print("Least likely substitution index is:", least_likely_subs)

# BLOSUM
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
#print(np.round(A_blos).astype(int))

Atot_blos = np.sum(np.sum(A_blos)).astype(float)

q = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    q[i][j] = A_blos[i][j] / Atot_blos

#print("\n")
#print(np.round(q, 2).astype(float))

# Calculation of relative frequency matrix  R  by  Ri,j  =  qi,j/πiπj

R_blos = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    R_blos[i][j] = q[i][j] / (pi[i] * pi[j])

#print(np.round(R_blos).astype(int))

#I converted zero values in R matrix to 1 to deal with -inf in S matrix
R2 = np.where(R_blos==0,1,R_blos)

#print(R_blos)  
#print(R2) 

# Calculation of log-odds matrix  S  for BLOSUM using  Si,j  = [10  log10   Ri,j]

S_blos = np.zeros((20,20))

for i in range(0,20):
  for j in range(0,20):
    S_blos[i][j] = 10 * np.log10(R2[i][j])

#print(S_blos)   
#print(np.round(S_blos).astype(int))

df_S_BLOS = pd.DataFrame(np.round(S_blos).astype(int), columns = ["A-0","R-1","N-2","D-3","C-4","Q-5","E-6","G-7","H-8","I-9","L-10","K-11","M-12","F-13","P-14","S-15","T-16","W-17","Y-18","V-19"], 
                         index = ["A-0","R-1","N-2","D-3","C-4","Q-5","E-6","G-7","H-8","I-9","L-10","K-11","M-12","F-13","P-14","S-15","T-16","W-17","Y-18","V-19"])
print("\n")
print("Log-odds matrix for BLOSUM:")
print("\n")
print(df_S_BLOS)
print("\n")

most_cons_PAM = np.where(S_blos == np.amax(S_blos[mask_diag]))
print("Most conserved aa index is:", most_cons_PAM)
least_cons_PAM = np.where(S_blos == np.amin(S_blos[mask_diag]))
print("Least conserved aa index is:", least_cons_PAM)

most_likely_subs = np.where(S_blos == np.amax(S_blos[mask_non_diag]))
print("Most likely substitution index is:", most_likely_subs)
least_likely_subs = np.where(S_blos == np.amin(S_blos[mask_non_diag]))
print("Least likely substitution index is:", least_likely_subs)
