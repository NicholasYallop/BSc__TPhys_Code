import time
import math
import numpy as np
import random as rand
runnum = 25
runtime = np.zeros([3, runnum+1])
absolute_error = 0 

def DET(A):
    DIM = A.shape[0]
    C = np.zeros((DIM,DIM))
    Determinant = 0
    if A.shape[0]>2 or A.shape[1]>2:
        for i in range(0,DIM):
            for j in range(0,DIM):
                irange = []
                jrange = []
                for k in range(0,DIM):
                    if k != i:
                        irange.append(k)
                    if k != j:
                        jrange.append(k)
                #print(jrange, irange)
                Aprime = A[irange,1:DIM]
                if j ==0:
                    Determinant += (DET(Aprime)[0])*A[i,0]*(-1)**(i)
                Amark =A[irange,:]
                Amark =Amark[:,jrange]
                C[j,i] = DET(Amark)[0]*(-1)**(i+j)
    else:
        Determinant = A[0,0]*A[1,1]-A[1,0]*A[0,1]
        #print(Determinant)
    if Determinant == 0:
        INV = np.dot(np.ones([dim,dim]), 999)
    else:
        INV = C/Determinant
    return(Determinant, C, INV)
     
for g in range(0,runnum):
    if g == 0:
        absolute_start = time.time()
    start_time = time.time()    
    "DEFINE INPUT"
    dim = 5
    Input=np.zeros((dim,dim))
    for i in range(0,dim):
        for j in range(0,dim):
             Input[i,j]=rand.randint(10,20)
             """
            if i==0 and j==0:
                Input[i,j]=-149
            if i==0 and j==1:
                Input[i,j]=-50
            if i==0 and j==2:
                Input[i,j]=-154
            if i==1 and j==0:
                Input[i,j]=537
            if i==1 and j==1:
                Input[i,j]=180
            if i==1 and j==2:
               Input[i,j]=546
            if i==2 and j==0:
                Input[i,j]=-27       
            if i==2 and j==1:
                Input[i,j]=-9
            if i==2 and j==2:
                Input[i,j]=-25
            """
    "Inverter"
    Answer = DET(Input)
    #print(Answer[2])
    error = 0
    errorMAT = np.dot(Answer[2],Input)
    for i in range(0,dim):
        for j in range(0,dim):
            if i == j:
                error += abs(errorMAT[i,j] - 1)
            else:
                error += abs(errorMAT[i,j])
    runtime[2, g] = time.time() - start_time
    if math.isnan(error):
        runtime[1, g] = 999
    else:
        runtime[1, g] = abs(error)
    runtime[0, g] = g
    absolute_error += runtime[1,g]
runtime[0, runnum] = -1
runtime[1, runnum] = absolute_error
runtime[2, runnum] = time.time() - absolute_start
#print(Answer)
print(runtime)
