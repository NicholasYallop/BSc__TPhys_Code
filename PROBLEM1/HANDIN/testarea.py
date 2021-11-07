import time
import math
import numpy as np
import random as rand

def eigen(A):
    dim=A.shape[0]
    Out = A
    RHold = np.identity(dim)
    n=0
    breaker = 0
    while n <= 10:
        R = np.identity(dim)
        RT = np.identity(dim)
        maxval = 0
        Imax = 0
        Jmax = 0
        for i in range(0,dim):
            for j in range(0,dim):
                if i != j :
                    if abs(Out[i,j])>abs(maxval):
                        maxval = Out[i,j]
                        Imax = i
                        Jmax = j
        #print(Imax,Jmax)
        if Imax == Jmax == 0:
            print("meme1")
            break
        if Out[Imax,Jmax] != 0:
            W = (Out[Jmax,Jmax]-Out[Imax,Imax])/(2*Out[Imax,Jmax])
        else:
            print("meme2")
            break
        if W>0:
            T = -W +math.sqrt(1+W**2)
        elif W<0:
            T = -W -math.sqrt(1+W**2)
        else:
            print("meme3", W)
            breaker += 1
            if breaker >= 25:
                break
            continue
        #print(W)
        #Out = rotate(Out, Imax, Jmax, T)
        #print(Out)
        C = (1+T*T)**(-0.5)
        S = C*T
        R[Imax,Imax] = C
        R[Jmax,Jmax] = C
        R[Imax,Jmax] = S
        R[Jmax,Imax] = -S
        RT[Imax,Imax] = C
        RT[Jmax,Jmax] = C
        RT[Imax,Jmax] = -S
        RT[Jmax,Imax] = S
        #print(S)
        #print(R)   
        Out = np.dot(RT,np.dot(Out,R))
        RHold = np.dot(RHold,R)
        if abs(S) < (10**-20):
            n+=1
        else:
            n=0
    eigenvects = RHold
    eigenvals = Out
    return(eigenvects, eigenvals)

def SVD(Input):
    start_time = time.time()               
    dimN=Input.shape[0] #defines the input matrix dimension
    dimM=Input.shape[1]
    "SVD approach"
    AAT = np.dot(Input,np.transpose(Input))
    ATA = np.dot(np.transpose(Input),Input)
    Ufull = eigen(AAT)
    Vfull = eigen(ATA)
    #Uorg = Organise(Ufull[0],Ufull[1])
    #Vorg = Organise(Vfull[0],Vfull[1])
    U = Ufull[0]
    Uprime = Ufull[1]
    V = Vfull[0] 
    Vprime = Vfull[1]      
    SigmaU = np.zeros([AAT.shape[0],ATA.shape[1]])
    SigmaV = np.zeros([AAT.shape[0],ATA.shape[1]])
    Sigma = np.zeros([AAT.shape[0],ATA.shape[1]])
    SigmaUT = np.zeros([ATA.shape[0],AAT.shape[1]])
    SigmaVT = np.zeros([ATA.shape[0],AAT.shape[1]])
    SigmaT = np.zeros([ATA.shape[0],AAT.shape[1]])
    mindim=min(AAT.shape[0],ATA.shape[1])
    Vtemp = V
    Utemp = U
    for k in range(0,dimM):
        TEST = np.dot(ATA, V[:,k])
        VecSum = 0
        for i in range(0,dimN):
            TEST[i]=TEST[i]/V[i,k]
            VecSum += TEST[i]
            if i != k:
                SigmaU[i,k]=Uprime[i,k]
        if Uprime[k,k] >= 0:
            SigmaU[k,k] = math.sqrt(Uprime[k,k])
        else:
            print(Uprime)
            break
        for j in range(0,dimN):
            if abs(Uprime[j,j] - VecSum/dimM) < 10**-3:
                U[:,j] = Utemp[:,k]
                
    for k in range(0,dimN):
        TEST = np.dot(AAT, U[:,k])
        VecSum = 0
        for i in range(0,dimN):
            TEST[i]=TEST[i]/U[i,k]
            VecSum += TEST[i]
            if i != k:
                SigmaV[i,k]=Vprime[i,k]
        if Vprime[k,k] >= 0:
            SigmaV[k,k] = math.sqrt(Vprime[k,k])
        else:
            print(SigmaV)
            break
        for j in range(0,dimN):
            if abs(Vprime[j,j] - VecSum/dimM) < 10**-3:
                V[:,j] = Vtemp[:,k]
        
    for i in range(0,mindim):
        SigmaUT[i,i]=1/SigmaU[i,i]
        SigmaVT[i,i]=1/SigmaV[i,i]
        for k in range(0,mindim):
            if i!= k:
                SigmaUT[k,i]=SigmaU[i,k]
                SigmaVT[k,i]=SigmaV[i,k]
    for i in range(0,mindim):
        for j in range(0,mindim):
            Sigma[i,j] = 0.5*(SigmaU[i,j]+SigmaV[i,j])
            SigmaT[i,j] = 0.5*(SigmaUT[i,j]+SigmaVT[i,j])
    fuggedup=0
    for i in range(0,dimN):
        if SigmaU[i,i] - SigmaV[i,i] > 0.01:
            fuggedup += 1
            #print(Uprime[i,i],Vprime[i,i])
    #if fuggedup != 0:
     #   print("AAT->ATA discrepency")
      #  print(Uprime, Vprime)
    VSTUT = np.dot(V,np.dot(SigmaT,np.transpose(U)))
    ErrorCheck = np.dot(VSTUT, Input)
    Error = 0
    for i in range(0,dimN):
        for j in range(0,dimN):
            if i!=j:
                Error += abs(ErrorCheck[i,j])
            else:
                Error += abs(ErrorCheck[i,j]-1)
    #print(Uprime,Vprime)
    return(U, V, Sigma, SigmaT, VSTUT, time.time()-start_time, Error)


print("Input:")
print(Input)
SV = SVD(Input)
print("U:")
print(SV[0])
print("V:")
print(SV[1])
print("Sigma:")
print(SV[2])
print("Inverse:")
print(SV[4])
print("Error:")
print(SV[6])