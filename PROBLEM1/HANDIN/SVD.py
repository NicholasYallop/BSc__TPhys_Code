import time
import math
import numpy as np
import random as rand


def multiply(A,B): #np.dot rework
    Output = np.zeros((A.shape[0],B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            if A.shape[1] == B.shape[0]:
                for k in range(0, A.shape[1]):
                    Output[i,j] += A[i,k]*B[k,j]
    return(Output)
    
def order(x,dim): #organizes singular values to size order
    SizeOrdered = 0 
    while SizeOrdered < dim-1:
        SizeOrdered = 0
        for i in range(0,dim-1):
            if abs(x[i])<abs(x[i+1]):
                Temp1=x[i]
                Temp2=x[i+1]
                x[i] = Temp2
                x[i+1] = Temp1
            else:
                SizeOrdered += 1

def rotate(A,i,j,t):
    dim = A.shape[0]
    B = np.zeros(A.shape)
    B[i,i]=(A[i,i] -2*t*A[i,j] +t*t*A[j,j])/(1+t*t)
    B[j,j]=(t*t*A[i,i] +2*t*A[i,j] +A[j,j])/(1+t*t)
    B[i,j]=((1-t*t)*A[i,j] + t*(A[i,i]-A[j,j]))/(1+t*t)
    B[j,i]=B[i,j]
    for k in range(0,dim):
        if k != i and k != j:
            B[i,k] = (A[i,k] - t*A[j,k])/(1+t*t)
            B[k,i] = B[i,k]
            B[j,k] = (t*A[i,k] + A[j,k])/(1+t*t)
            B[k,j] = B[k,j]
            for l in range(0,dim):
                if l != i and l != j:
                    B[k,l] = A[k,l]
    return(B)
    
"Jacobi Eigenfinder"
def eigen(A):
    dim=A.shape[0]
    Out = A
    RHold = np.identity(dim)
    n=0
    breaker = 0
    while n != 2:
        R = np.identity(dim)
        RT = np.identity(dim)
        maxval = 0
        Imax = 0
        Jmax = 0
        Iblack, Jblack = 0,0
        for i in range(0,dim):
            for j in range(0,dim):
                if i != j :
                    if i != Iblack or j != Jblack:
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
            Iblack = Imax
            Jblack = Jmax
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
        #print(R)   
        Out = np.dot(RT,np.dot(Out,R))
        RHold = np.dot(RHold,R)
        if abs(S) < (10**-18):
            n+=1
        else:
            n=0
    eigenvects = RHold
    eigenvals = Out
    return([eigenvects, eigenvals])
  
def Organise(M, Mprime): #aligns eigenvectors to correct singular value
    dim = M.shape[0]
    Morganised = 0
    MEigs = np.zeros(dim)
    for i in range(0,dim):
        MEigs[i]=Mprime[i,i]
    while Morganised < dim-1:
        Morganised = 0
        for i in range(0,dim-1):
            if Mprime[i,i] > Mprime[i+1,i+1]:
                Mtemp = np.copy(Mprime[:, i])
                Mprime[:, i] = Mprime[:, i+1]
                Mprime[:, i+1] = Mtemp
                Hold = np.zeros((dim,2))
                for k in range(0,dim):
                    if k == dim-1:
                        Hold[k,0] = Mprime[0,i]
                    else:
                        Hold[k,0] = Mprime[k+1,i]
                    if k == 0:
                        Hold[k,1] = Mprime[dim-1,i+1]
                    else:
                        Hold[k,1] = Mprime[k-1,i+1]
                Mprime[:, i] = Hold[:,0]
                Mprime[:, i+1] = Hold[:,1]
            else:
                Morganised += 1
    MtempMAT = M
    M=np.zeros((dim,dim))
    for i in range(0,dim):
        for j in range(0,dim):
            if Mprime[i,i] == MEigs[j]:
                M[:,j] = MtempMAT[:,i]
    return(M, Mprime)                 

runnum = 1
runtime = np.zeros([3, runnum+1])
absolute_error = 0                 
for g in range(0,runnum):
    if g==0:
        absolute_start = time.time()
    start_time = time.time()               
    dimN=3 #defines the input matrix dimension
    dimM=3
    
    "DEFINE INPUT"
    Input=np.zeros((dimM,dimN))
    for i in range(0,dimM):
        for j in range(0,dimN):
             Input[i,j]=rand.randint(1,10)
             #Input[i,j]=rand.uniform(-1,1)
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
    
    "SVD approach"
    AAT = np.dot(Input,np.transpose(Input))
    ATA = np.dot(np.transpose(Input),Input)
    U, Uprime = eigen(AAT)
    V, Vprime = eigen(ATA)
    U, Uprime = Organise(U, Uprime)
    V, Vprime = Organise(V, Vprime)        
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
    if fuggedup != 0:
        print("ya dun fugged up")
        
    USVT = np.dot(U,np.dot(Sigma,np.transpose(V)))
    VSTUT = np.dot(V,np.dot(SigmaT,np.transpose(U)))
    UUpUT = np.dot(U,np.dot(Uprime,np.transpose(U)))
    VVpVT = np.dot(V,np.dot(Vprime,np.transpose(V)))
    #print(Uprime)
    TEST = np.dot(ATA, V[:,0]) 
    for i in range(0,dimN):
        TEST[i]=TEST[i]/V[i,0]
    ErrorCheck = np.dot(VSTUT, Input)
    Error = 0
    for i in range(0,dimN):
        for j in range(0,dimN):
            if i!=j:
                Error += abs(ErrorCheck[i,j])
            else:
                Error += abs(ErrorCheck[i,j]-1)
    runtime[2, g] = time.time() - start_time
    if math.isnan(Error):
        runtime[1, g] = 999
    else:
        runtime[1, g] = abs(Error)
    runtime[0, g] = g
    absolute_error += runtime[1,g]
    print(ErrorCheck)
runtime[0, runnum] = -1
runtime[1, runnum] = absolute_error
runtime[2, runnum] = time.time() - absolute_start
print(runtime)
#print(UUpUT, AAT)
#print(VVpVT, ATA)
#print(USVT, Input)
#print(VSTUT, np.linalg.inv(Input))
#print(np.dot(np.linalg.inv(Input),Input))