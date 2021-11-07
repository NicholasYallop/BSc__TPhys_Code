import cmath
import math
import numpy as np
import random as rand
def multiply(A,B):
    Output = np.zeros((A.shape[0],B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            if A.shape[1] == B.shape[0]:
                for k in range(0, A.shape[1]):
                    Output[i,j] += A[i,k]*B[k,j]
    return(Output)
dim=3
iterations = 100
attempts = 100
Input=np.ones((dim,dim))
"define input for tests:"
for i in range(0,dim):
    for j in range(0,dim):
        Input[i,j]= float(rand.randint(0,10))
"Jacobi Eigenfinder"
A = Input
RHold = np.identity(dim)
for n in range(0,100):
    R = np.identity(dim)
    RT = np.identity(dim)
    maxval = 0
    Imax = 0
    Jmax = 0
    for i in range(0,dim):
        for j in range(0,dim):
            if i!= j:
                if abs(A[i,j])>maxval:
                    maxval = A[i,j]
                    Imax = i
                    Jmax = j
    if Imax == Jmax == 0:
        print("meme1")
        break
    if A[Imax,Jmax] != 0:
        W = (A[Jmax,Jmax]-A[Imax,Imax])/(2*A[Imax,Jmax])
    else:
        print("meme2")
        break
    if W>0:
        T = -W +math.sqrt(1+W**2)
    elif W<0:
        T = -W -math.sqrt(1+W**2)
    else:
        print("meme3")
        continue
    C = (1+T*T)**(-0.5)
    S = C*T
    R[Imax,Imax] = C
    R[Jmax,Jmax] = C
    R[Imax,Jmax] = S
    R[Jmax,Imax] = -S
    RT[Jmax,Jmax] = C
    RT[Imax,Imax] = C
    RT[Jmax,Imax] = S
    RT[Imax,Jmax] = -S
    RHold = multiply(RHold, R)
    #print(RHold)
    Outtemp = multiply(RT,A)
    Out = multiply(Outtemp, R)
    A = Out
print(np.linalg.eig(Input))
print(multiply(multiply(RHold,A),np.transpose(RHold)))
                
    
    