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
    Holdingon = 0
    while n <= 3:
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
            #print("meme1")
            break
        if Out[Imax,Jmax] != 0:
            W = (Out[Jmax,Jmax]-Out[Imax,Imax])/(2*Out[Imax,Jmax])
        else:
            #print("meme2")
            break
        if W>0:
            T = -W +math.sqrt(1+W**2)
        elif W<0:
            T = -W -math.sqrt(1+W**2)
        else:
            #print("meme3", W)
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
        if abs(S) < (10**-20) or Holdingon > 200:
            n+=1
        else:
            n=0
            Holdingon += 1
    eigenvects = RHold
    eigenvals = Out
    return(eigenvects, eigenvals)
  
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

def DET(A):
    start_time = time.time()
    dim = A.shape[0]
    C = np.zeros((dim,dim))
    Determinant = 0
    if A.shape[0]>2 or A.shape[1]>2:
        for i in range(0,dim):
            for j in range(0,dim):
                irange = []
                jrange = []
                for k in range(0,dim):
                    if k != i:
                        irange.append(k)
                    if k != j:
                        jrange.append(k)
                #print(jrange, irange)
                Aprime = A[irange,1:dim]
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
    error = 0
    errorMAT = np.dot(INV,A)
    for i in range(0,dim):
        for j in range(0,dim):
            if i == j:
                error += abs(errorMAT[i,j] - 1)
            else:
                error += abs(errorMAT[i,j])
    return(Determinant, C, INV, time.time()-start_time, error)
  
def LU(Input):
    start_time = time.time()
    dim=Input.shape[0] #defines the input matrix dimension
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
                if U[i,i] == 0:
                    L[k,i] = 99999
                else:
                    L[k,i]=(Input[k,i]-Ldiff)/U[i,i]  #Solves for L cell
    "LU Checker"
    LUCheck=0
    for i in range(0,dim):
        for j in range(0,dim):
            LUCheck += np.dot(L,U)[i,j] - Input[i,j]
    #if LUCheck >= 10**(-10):
        #print("LU reconstruction inaccurate", LUCheck)
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
            if E[i,i] == 0:
                BU[i,:] = np.ones(dim)*99999
            else:
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
            if E[i,i] == 0:
                BU[0,i]=99999 
            else:
                BU[0,i]=BU[0,i]/E[i,i]
        UI[k,k:dim] = BU[0,:]
    Inverse = np.dot(UI,LI)
    Ident = np.dot(Inverse, Input)
    Error = 0
    for i in range(0,dim):
        for j in range(0,dim):
            if i!=j:
                Error += abs(Ident[i,j])
            else:
                Error += abs(Ident[i,j]-1)
    return(L, U, LI, UI, np.dot(UI, LI), time.time()-start_time, Error)

def SVD(Input):
    start_time = time.time()               
    dimN=Input.shape[0] #defines the input matrix dimension
    dimM=Input.shape[1]
    "SVD approach"
    AAT = np.dot(Input,np.transpose(Input))
    ATA = np.dot(np.transpose(Input),Input)
    Ufull = eigen(AAT)
    Vfull = eigen(ATA)
    Uorg = Organise(Ufull[0],Ufull[1])
    Vorg = Organise(Vfull[0],Vfull[1])
    U = Uorg[0]
    Uprime = Uorg[1]
    V = Vorg[0] 
    Vprime = Vorg[1]      
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
            #print(Uprime)
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
            #print(SigmaV)
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

def INV(Input):
    La = LU(Input)
    if La[6]<10**-10:
        #print(La[6], "LU")
        return(La[4])
    else:
        SV = SVD(Input)
        if SV[6]<10**-5:
            #print(SV[6], "SVD")
            if La[6] < SV[6]:
                return(La[4])
            else:
                return(SV[4])
        else:
            De = DET(Input)
            #print(De[4], DET)
            return(De[2])
    
def TENSION(position):
    gravity = np.array([0,0,-9.81])
    mass = 70
    drum1 = np.array([0,0,8])
    drum2 = np.array([0,15,8])
    drum3 = np.array([8,7.5,8])
    theta1 = drum1 - position
    theta2 = drum2 - position
    theta3 = drum3 - position
    theta1 = theta1/math.sqrt(np.dot(theta1,theta1))
    theta2 = theta2/math.sqrt(np.dot(theta2,theta2))
    theta3 = theta3/math.sqrt(np.dot(theta3,theta3))
    theta = np.zeros([3,3])
    theta[:,0] = theta1
    theta[:,1] = theta2
    theta[:,2] = theta3
    weight = np.dot(mass, gravity)
    return(np.dot(INV(theta), -weight))

def MAXTEN(xsample,ysample,zsample,xhi,yhi,zhi, drumnumber):
    xlow, ylow, zlow = 0,0,0
    imaxten = 0
    jmaxten, jmaxtenlast = 0, 0
    kmaxten, kmaxtenlast = 0, 0
    location = np.zeros(3)
    location1 = location
    location2 = location1
    for i in range(0,xsample):
        for j in range(0,ysample):
            for k in range(0,zsample):
                maxten = 0
                sample1 = [xlow+(i/xsample)*xhi,ylow+(j/ysample)*yhi,zlow+(k/zsample)*zhi]
                sample2 = [xlow+((i+1)/xsample)*xhi,ylow+((j+1)/ysample)*yhi,zlow+((k+1)/zsample)*yhi]
                TEN1 = TENSION(sample1)
                TEN2 = TENSION(sample2)
                for i in range(0,3):
                    if TEN1[i]<0 or TEN2[i]<0:
                        TEN1, TEN2 = np.zeros(3), np.zeros(3)
                print(TEN1)
                if  TEN1[drumnumber] >= TEN2[drumnumber]:
                    break
                else:
                    if (TEN2[drumnumber]) > (maxten):
                        maxten = TEN2[drumnumber]
                if (maxten) > (kmaxten):
                    kmaxten = maxten
                    location = sample2
                print(i,j,k)
            if (kmaxten) < (kmaxtenlast):
                break
            else:
                kmaxtenlast = kmaxten
            if (kmaxten) > (jmaxten):
                jmaxten = kmaxten
                location1 = location
        if (jmaxten) < (jmaxtenlast):
            break
        else:
            jmaxtenlast = jmaxten
        if (jmaxten) > (imaxten):
            imaxten = jmaxten
            location2 = location1
    return(imaxten, location2)
    
print("Nick Yallop - ny16281")
print("Excercise 1 - matrix inversion and wire tension")
finish = 0
while finish != 1:
    print("[1] - Cramer's Rule")
    print("[2] - LU Decomposition")
    print("[3] - Singular Value Decomposition")
    print("[4] - Assigned Physical Problem")
    METHcorrectinput = 0
    while METHcorrectinput == 0:
        METHselection = input("Method Selection:")
        if METHselection == "1" or METHselection == "2" or METHselection == "3" or METHselection == "4":
            METHcorrectinput=1
        else:
            print("come on now, pick a number from 1 to 4")
    if METHselection != "4":
        print("DEFINE INPUT MATRIX")
        dim = input("Number of Dimensions:")
        DIMcorrect = 0
        while DIMcorrect == 0:
            try:
                DIMinput = int(dim)
            except ValueError:
                print("Pick an integer")
                dim = input("Number of Dimensions:")
            else:
                DIMcorrect = 1
        dim = int(dim)
        Input = np.zeros([dim,dim])
        print("[1] - Random Cells")
        print("[2] - Random, Near-Singular")
        print("[3] - User Input")
        INPcorrectinput = 0
        while INPcorrectinput == 0:
            INPselection = input("Selection:")
            if INPselection == "1" or INPselection == "2":
                INPcorrectinput=1
            else:
                print("come on now, pick 1, 2 or 3")
        
        if INPselection == "1" or INPselection =="2":
            print("Uniform Random or Integer Random?")
            print("[1] - Uniform")
            print("[2] - Integer")
            RANcorrectinput = 0
            while RANcorrectinput == 0:
                RANselection = input("Selection:")
                if RANselection == "1" or RANselection == "2":
                    RANcorrectinput=1
                else:
                    print("come on now, pick a number from 1 to 2")
            print("BOUNDS FOR RANDOM CELLS")
            lower = input("Lower Bound:")
            LOWcorrect = 0
            while LOWcorrect == 0:
                try:
                    LOWinput = float(lower)
                except ValueError:
                    print("Pick a number")
                    lower = input("Lower Bound:")
                else:
                    LOWcorrect = 1
            upper = input("Upper Bound:")
            UPPcorrect = 0
            while UPPcorrect == 0:
                try:
                    UPPinput = float(upper)
                except ValueError:
                    print("Pick a number")
                    upper = input("Upper Bound:")
                else:
                    UPPcorrect = 1
            lower = float(lower)
            upper = float(upper)
            for i in range(0,dim):
                for j in range(0,dim):
                    if RANselection == "1":
                        Input[i,j] = rand.uniform(lower,upper)
                    if RANselection == "2":
                        Input[i,j] = rand.randint(lower,upper)
        if INPselection == "2":
            Input[:,dim-1]=np.zeros(dim)
            for i in range(0,dim-1):
                Input[:,dim-1] += Input[:,i]
            Input[:,dim-1] = Input[:,dim-1]/math.sqrt(np.dot(Input[:,dim-1],Input[:,dim-1]))
            #for i in range(0,dim):
             #   Input[i,dim-1] += rand.uniform(0.000001*lower,0.000001*upper) 
        if INPselection == "3":
            for i in range(0,dim):
                for j in range(0,dim):
                    print("cell [",i,",",j,"]")
                    InputHold = input("Value:")
                    UINcorrect = 0
                    while UINcorrect == 0:
                        try:
                            UINinput = float(InputHold)
                        except ValueError:
                            print("Pick a number")
                            print("cell [",i,",",j,"]")
                            InputHold = input("Value:")
                        else:
                            UINcorrect = 1
                    Input[i,j] = InputHold
        
        if METHselection == "1":
            output = DET(Input)
            print("Input:")
            print(Input)
            print("Determinant:", output[0])
            print("Inverse:")
            print(output[2])
            print("Time:", output[3])
            print("Error:", output[4])
        if METHselection == "2":
            output = LU(Input)
            print("Input:")
            print(Input)
            print("L:")
            print(output[0])
            print("U:")
            print(output[1])
            print("Inverse:")
            print(output[4])
            print("Time:", output[5])
            print("Error:", output[6])
        if METHselection == "3":
            output = SVD(Input)
            print("Input:")
            print(Input)
            print("U:")
            print(output[0])
            print("V:")
            print(output[1])
            print("Sigma:")
            print(output[2])
            print("Inverse:")
            print(output[4])
            print("Time:", output[5])
            print("Error:", output[6])
    if METHselection == "4":
        print("[1] - Tension Finder")
        print("[2] - Max Tension Location")
        PHYcorrectinput = 0 
        while PHYcorrectinput == 0:
            PHYselection = input("Selection:")
            if PHYselection == "1" or PHYselection == "2":
                PHYcorrectinput=1
            else:
                print("come on now, pick 1 or 2")
        if PHYselection == "1":
            xpos = input("Distance from front of stage (m):")
            LOCxcorrect = 0
            while LOCxcorrect == 0:
                try:
                    LOCxinput = float(xpos)
                except ValueError:
                    print("Pick a number")
                    xpos = input("Distance from front of stage (m):")
                else:
                    if float(xpos) <= 8:
                        LOCxcorrect = 1
                    else:
                        print("outside bounds")
                        xpos = "apple"
            xpos = float(xpos)
            ypos = input("Distance from left of stage (m):")
            LOCycorrect = 0
            while LOCycorrect == 0:
                try:
                    LOCyinput = float(ypos)
                except ValueError:
                    print("Pick a number")
                    ypos = input("Distance from left of stage (m):")
                else:
                    if float(ypos) <= 15:
                        LOCycorrect = 1
                    else:
                        print("outside bounds")
                        ypos = "apple"
            ypos = float(ypos)
            zpos = input("Distance above ground level (m):")
            LOCzcorrect = 0
            while LOCzcorrect == 0:
                try:
                    LOCzinput = float(zpos)
                except ValueError:
                    print("Pick a number")
                    zpos = input("Distance above ground level (m):")
                else:
                    if float(zpos) <= 8:
                        LOCzcorrect = 1
                    else:
                        print("outside bounds")
                        zpos = "apple"
            zpos = float(zpos)
            output = TENSION([xpos, ypos, zpos])
            print("Tension in String 1:", output[0])
            print("Tension in String 2:", output[1])
            print("Tension in String 3:", output[2])
        if PHYselection == "2":
            xsamp = input("X sample number:")
            SAMPxcorrect = 0
            while SAMPxcorrect == 0:
                try:
                    SAMPxinput = int(xsamp)
                except ValueError:
                    print("Pick a integer")
                    xsamp = input("X sample number:")
                else:
                    SAMPxcorrect = 1
            xsamp = int(xsamp)
            ysamp = input("Y sample number:")
            SAMPycorrect = 0
            while SAMPycorrect == 0:
                try:
                    SAMPyinput = int(ysamp)
                except ValueError:
                    print("Pick a integer")
                    ysamp = input("Y sample number:")
                else:
                    SAMPycorrect = 1
            ysamp = int(ysamp)
            zsamp = input("Z sample number:")
            SAMPzcorrect = 0
            while SAMPzcorrect == 0:
                try:
                    SAMPzinput = int(zsamp)
                except ValueError:
                    print("Pick a integer")
                    zsamp = input("Z sample number:")
                else:
                    SAMPzcorrect = 1
            zsamp = int(zsamp)
            print("Select Drum - 1, 2 or 3")
            DRUMcorrectinput = 0 
            while DRUMcorrectinput == 0:
                DRUMselection = input("Selection:")
                if DRUMselection == "1" or DRUMselection == "2" or DRUMselection == "3":
                    DRUMcorrectinput=1
                else:
                    print("come on now, pick 1, 2 or 3")
            DRUMno = int(DRUMselection)-1
            print(MAXTEN(xsamp,ysamp,zsamp,8,15,8,DRUMno))
    print("Finished?")
    print("[1] - Yes, quit")
    print("[2] - No, from the top")            
    QUITcorrectinput = 0 
    while QUITcorrectinput == 0:
        QUITselection = input("Selection:")
        if QUITselection == "1" or QUITselection == "2":            
            QUITcorrectinput=1
        else:
            print("come on now, pick 1 ro 2")
    if QUITselection == "1":
        finish = 1