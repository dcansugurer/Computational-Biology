# Import required libraries

import pandas as pd
import numpy as np

font = {'family'  :   'Dejavu Sans',
        'weight'  :   'normal',
        'size'    :   20}

matplotlib.rc('font', **font)

# SQ1 = human insulin amino acid sequence
# SQ2 = mouse insulin amino acid sequence

SQ1 = 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN'
SQ2 = 'MALLVHFLPLLALLALWEPKPTQAFVKQHLCGPHLVEALYLVCGERGFFYTPKSRREVEDPQVEQLELGGSPGDLQTLALEVARQKRGIVDQCCTSICSLYQLENYCN'

L1 = len(SQ1)
L2 = len(SQ2)

alpha = 2

beta = -1

omega = -2

S = np.zeros((L1+1,L2+1),'int')
D = list()
D.append([])

for i in range(0,L1):
    D.append([])
    D[i+1].append([])
    for j in range(0,L2):
        if SQ1[i]==SQ2[j]:
            S[i+1][j+1] = alpha
        else:
            S[i+1][j+1] = beta
        D[i+1].append([])
        
A = np.zeros((L1+1,L2+1),'int')

for i in range(1,L1+1):
    for j in range(1,L2+1):
        score_from_above = A[i-1][j] + omega
        score_from_left = A[i][j-1] + omega
        score_from_diag = A[i-1][j-1] + S[i][j]
        
        best_score = np.max([score_from_above, score_from_left, score_from_diag])
        A[i][j] = best_score
        if score_from_diag == best_score:
            D[i][j].append('d')
        if score_from_above == best_score:
            D[i][j].append('a')
        if score_from_left == best_score:
            D[i][j].append('l')
            
print(A)

ind_max_lastrow = np.argmax(A[L1,:])
ind_max_lastcol = np.argmax(A[:,L2])

if A[L1,ind_max_lastrow] >= A[ind_max_lastcol,L2]:
    i_end = L1
    j_end = ind_max_lastrow
else:
    i_end = ind_max_lastcol
    j_end = L2


ci = i_end - 1
cj = j_end - 1

SQ1_aln = SQ1[ci: L1]
SQ2_aln = SQ2[cj: L2]

M_aln = '|' if SQ1[ci] == SQ2[cj] else ' '

while (ci > 0 and cj > 0):

    if D[ci][cj][0] == 'd':
        ci -= 1
        cj -= 1
        SQ1_aln = SQ1[ci] + SQ1_aln
        SQ2_aln = SQ2[cj] + SQ2_aln

    else:
        if D[ci][cj][0] == 'a':
            ci -= 1
            SQ1_aln = SQ1[ci] + SQ1_aln
            SQ2_aln = '_' + SQ2_aln
        else:
            cj -= 1
            SQ1_aln = '_' + SQ1_aln
            SQ2_aln = SQ2[cj] + SQ2_aln

    M_aln = '|' + M_aln if SQ1[ci] == SQ2[cj] else ' ' + M_aln

print(SQ1_aln)
print(M_aln)
print(SQ2_aln)
