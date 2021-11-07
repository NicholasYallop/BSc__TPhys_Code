import time
import math
import numpy as np
import random as rand
runnum = 25
runtime = np.zeros([3, runnum+1])
absolute_error = 0
for g in range(0,runnum):
    if g == 0:
        absolute_start = time.time()
    start_time = time.time()
    dim=5 #defines the input matrix dimension
    Input=np.zeros((dim,dim))
    "define random input for tests:"
    for i in range(0,dim):
        for j in range(0,dim):
            Input[i,j]=rand.randint(1,10)
    #print(Input)
    "LU approach:"
    L=np.zeros((dim,dim))#generate empty L
    U=np.zeros((dim,dim))#generate empty U
    Ldiff = np.zeros(dim)
    Udiff = np.zeros(dim)
    for i in range(0,dim): #runthrough iteration
        L[i,i] = 1 #L diagonal = 1s
        for k in range(i,dim): #parses cells from a row below/column in to the limit of the matrix
            Udiff = 0
            Ldiff = 0
            for m in range(0,i):
                Udiff += L[i,m]*U[m,k] #Difference between input and U value
                Ldiff += L[k,m]*U[m,i] #Difference between input and L denominator
            if i > k:
                U[i,k] = 0
                L[k,i] = 0               
            else:
                U[i,k]=Input[i,k]-Udiff #Solves for U cell
                L[k,i]=(Input[k,i]-Ldiff)/U[i,i]  #Solves for L cell
    "LU Checker"
    LUCheck=0
    for i in range(0,dim):
        for j in range(0,dim):
            LUCheck += np.dot(L,U)[i,j] - Input[i,j]
    if LUCheck < 10**(-10):
        print("success")
    else:
        print(LUCheck)
    "Inverter"
    AL = L - np.identity(dim)
    LI = np.identity(dim)
    Lpowers = np.zeros((dim, dim, dim))
    Lpowers[0,:,:] = AL
    for i in range(1,dim):
        Lpowers[i,0:dim,0:dim] = np.dot(Lpowers[i-1,:,:],AL)
    for i in range(0,dim):
        LI += Lpowers[i,:,:]*((-1)**(i+1))
    UI = np.zeros(U.shape)
    for k in range(0,dim):
        E = U[k:dim,k:dim]
        BU = np.zeros((E.shape))
        for i in range(0,dim-k):
            BU[i,:] = E[i,:]/E[i,i]
        BU = BU - np.identity(dim-k)
        BU[0,0] = -1
        dummy = 1
        runcounter = 0
        while dummy == 1 and runcounter < 100:
            MatSum = 0
            BU = np.dot(BU,BU)
            for i in range(1,dim-k):
                for j in range(0,dim-k):
                    MatSum += abs(BU[i,j])
            if abs(MatSum) < 10**(-20):
                dummy=0
            else:
                runcounter += 1
        for i in range(0,dim-k):
            BU[0,i]=BU[0,i]/E[i,i]
        UI[k,k:dim] = BU[0,:]
    Inverse = np.dot(UI,LI)
    Ident = np.dot(Inverse, Input)
    Error = 0
    for i in range(0,dim):
        for j in range(0,dim):
            if i!=j:
                Error += Ident[i,j]
            else:
                Error += Ident[i,j]-1
    runtime[2, g] = time.time() - start_time
    if math.isnan(Error):
        runtime[1, g] = 0
    else:
        runtime[1, g] = abs(Error)
    runtime[0, g] = g
    absolute_error += runtime[1,g]
runtime[0, runnum] = -1
runtime[1, runnum] = absolute_error
runtime[2, runnum] = time.time() - absolute_start
#print(Error)

print(runtime)
    

    
            
