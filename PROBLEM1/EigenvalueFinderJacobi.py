import cmath
import math
import numpy as np
import random as rand
dim=3
iterations = 100
attempts = 100
Input=np.zeros((dim,dim))
"define random input for tests:"
for i in range(0,dim):
    for j in range(0,dim):
        Input[i,j]=rand.randint(0,10)
print(np.linalg.eig(Input)[0])
"eigenfinder"
eigenvects = np.zeros((dim,dim))
x = np.ones((dim, iterations))
A = Input
Output = np.zeros((dim,dim))
V = np.identity(dim)
XP=np.zeros(dim)
XQ=np.zeros(dim)
IBlacklist = [0]
JBlacklist = [0]
for n in range(0,100):
    R = np.identity(dim)
    RT = np.identity(dim)
    location = [0,0]
    maxval = 0
    for i in range(0,dim):
        for j in range (0,dim):
            if i != j:
                if i != IBlacklist and j != JBlacklist:
                    #print(abs(A[i,j]))
                    if abs(A[i,j])>maxval:
                        maxval = A[i,j]
                        I=i
                        J=j
    print(maxval)
    location=[I,J]
    if maxval == 0:
        break
        #continue
    if A[I,J] != 0:
        w = (A[J,J]-A[I,I])/(2*A[I,J])
    else:
        break
    if w>0:
        t = -w +math.sqrt(1+w*w)
    if w<0:
        t = -w -math.sqrt(1+w*w)
    Output=A
    #print(1/math.sqrt(1+t*t))
    Output[I,J] = 0
    Output[J,I] = 0
    Output[I,I] = ((A[I,I] + t*A[J,J]) - 2*t*A[I,J])/(1+t*t)
    Output[J,J] = ((A[I,I] + t*A[J,J]) + 2*t*A[I,J])/(1+t*t)
    for j in range(0,dim):
        if j != I and j!= J:
            Output[j,I] = (A[j,I] - t*A[j,J])/math.sqrt(1+t*t)
            Output[I,j] = Output[j,I]
            Output[j,J] = (A[j,J] + t*A[j,I])/math.sqrt(1+t*t)
            Output[J,j] = Output[j,J]        
    for j in range(0,dim):
        XP[j] = V[j,I]
        XQ[j] = V[j,J]
    for j in range(0,dim):
        V[j,I] = (XP[j]-t*XQ[j])/math.sqrt(1+t*t)
        V[j,J] = (t*XP[j]+XQ[j])/math.sqrt(1+t*t)
    breakcounter=0
    for i in range(0,dim):
        for j in range(0,dim):
            if i != j:
                if A[i,j] == 0:
                    breakcounter+=1
    if breakcounter == dim**2 - dim:
        break
        print("2")
    IBlacklist.append(I)
    JBlacklist.append(J)
    A = Output
    print(Output)
print(maxval)
#print(d)