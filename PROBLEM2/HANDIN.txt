import math as m 
import numpy as np 
import random as ran
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def V(x, y, h, ChargeDensity):
    return((Vmat[x-1,y] + Vmat[x+1,y] + Vmat[x,y+1] + Vmat[x,y-1] + h*h*ChargeDensity[x,y])/4)
    
def Vcheck(x, y, h, ChargeDensity):
    return((Vmat[x-1,y] + Vmat[x+1,y] + Vmat[x,y-1] + Vmat[x,y+1] - 4*Vmat[x,y])/(h*h) + ChargeDensity[x,y])
    
def RotMat(alpha, k, h, dim):
    a = - alpha/(h*h)
    b = (1/k) + (2*alpha)/(h*h)
    bp = (1/k) + (alpha)/(h*h)
    c = a
    RotMat = np.zeros((dim,dim))
    for i in range(0,dim):
        RotMat[i,i] = b                
        if i != 0:
            RotMat[i-1,i] = a
        if i == 0:
            RotMat[i+1,i] = c
            RotMat[i,i] = bp
        if i != dim-1:
            RotMat[i+1,i] = c
        if i == dim-1:
            RotMat[i,i] = bp
            RotMat[i-1,i] = a
    return(RotMat)
    
def ABS(vect):
    ABS = 0
    for i in range(0,vect.shape[0]):
        ABS += abs(vect[i])
    return(ABS)
        

def TimeProgress(TotalTime, Rod, RotMat, alphaair, alphaice, k, Tsteps, ice=0, oven=0):
    for t in range(1,Tsteps):
        RotMatInv = np.linalg.inv(RotMat)
        d = (1/k)*Rod    
        Rod = np.dot(RotMatInv, d)
        for i in range(0,oven):
            Rod[i] += (1000-Rod[i])*k*alphaair
        for i in range(0,ice):
            Rod[dim-i-1] += (-Rod[dim-i-1])*k*alphaice
        ABS = 0
        for i in range(0,dim):
            ABS += Rod[i]
        ABSPrint[t]=ABS
        #print(ABS)
        RodPrint[:,t] = Rod[0:dim]
    """
    for t in range(1,Tsteps):
        Rod = StepUp(Rod, RotMat, alphaair, alphaice, k, ice=0, oven=0)[0]
        RodPrint[:,t] = Rod
        ABS = 0
        for i in range(0,dim):
            ABS += abs(Rod[i])
        ABSPrint[t] = ABS
    """
    return(RodPrint, ABSPrint)
        
def Relax(ChargeDensity, Vmat, h, dim):
    VmatNew = Vmat
    n = 0
    Break = 0
    while Break == 0:
        Break = 1
        n += 1
        #print(n)
        for i in range(1,dim-1):
            for j in range(1,dim-1):
                VmatNew[i,j] = V(i, j, h, ChargeDensity)
                if Vmat[i,j] != 0:
                    if abs((Vmat[i,j] - VmatNew[i,j])/Vmat[i,j]) > 10**-7:
                        Break = 0
                else:
                    if abs(Vmat[i,j] - VmatNew[i,j]) > 10**-8:
                        Break = 0
        Vmat = VmatNew
        #print(Vmat[int(dim/2),:])
        #Vmat = 0.7*VmatNew + 0.3*Vmat
        if Break == 1:
            for i in range(1,dim-1):
                for j in range(1,dim-1):
                    HOLD = Vcheck(i,j,h,ChargeDensity)
                    if HOLD > 0.0001:
                        Break = 0
                    else:
                        Vcheckmat[i-1,j-1] = abs(HOLD)
        if n > 10000:
            Break = 1
    return(Vmat)
        
def Display(dim, T, k, RodPrint, inv=0):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Make data.
    X = np.arange(0, dim, 1)
    Y = np.arange(0, T, k)
    if inv == 0:
        X, Y = np.meshgrid(X, Y)
        Z = np.transpose(RodPrint)
    else:
        X, Y = np.meshgrid(Y, X)
        Z = RodPrint
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    # Customize the z axis.
    #ax.set_zlim(0, 100)
    ax.set_ylabel("Time (s)")
    if inv == 0:
        ax.set_xlabel("Distance Along Rod (m)")
        ax.set_zlabel("Temperature (deg C)")
        ax.set_ylabel("Time (s)")
    if inv == 1:
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Distance Along Rod (m)")
        ax.set_zlabel("Temperature (deg C)")
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

                        

alpha = 59/(7900*450)
alphaair = 0.2435/(1006*1.276)
alphaice = 2.22/(916.2*2.05)
              


print("Nick Yallop - ny16281")
print("Excercise 2 - Partial Differential Equations")
finish = 0
while finish != 1:
    print("[1] - Laplace Solver")
    print("[2] - Poisson Parallel Plate Solver")
    print("[3] - Heat Flow")
    METHcorrectinput = 0
    while METHcorrectinput == 0:
        METHselection = input("Method Selection:")
        if METHselection == "1" or METHselection == "2" or METHselection == "3":
            METHcorrectinput=1
        else:
            print("come on now, pick a number from 1 to 3")
    if METHselection != "3":
        print("DEFINE CONDITIONS")
        dim = input("Length of Matrix Dimensions:")
        DIMcorrect = 0
        while DIMcorrect == 0:
            try:
                DIMinput = int(dim)
            except ValueError:
                print("Pick an integer")
                dim = input("Length of Matrix Dimensions:")
            else:
                DIMcorrect = 1
        dim = int(dim)
        h = input("Length of Site (m):")
        Hcorrect = 0
        while Hcorrect == 0:
            try:
                Hinput = float(h)
            except ValueError:
                print("Pick an integer")
                h = input("Length of Site (m):")
            else:
                Hcorrect = 1
        h = float(h)
        Vmat = np.zeros((dim, dim))
        VmatNew = np.zeros((dim, dim))
        Vcheckmat = np.zeros((dim-2, dim-2))
        ChargeDensity = np.zeros((dim, dim))
        Input = np.zeros([dim,dim])
        print("Boundary Values")
        print("[1] - Random Values")
        print("[2] - Constant Value")
        print("[3] - i/(j+1) (normalised)")
        print("[4] - sqrt(i^2 + j^2) (normalised)")
        INPcorrectinput = 0
        while INPcorrectinput == 0:
            INPselection = input("Selection:")
            if INPselection == "1" or INPselection == "2" or INPselection == "3" or INPselection == "4":
                INPcorrectinput=1
            else:
                print("come on now, pick 1, 2, 3 or 4")
        if INPselection == "1":
            print("Uniform Random or Integer Random?")
            print("[1] - Uniform")
            print("[2] - Integer")
            RANcorrectinput = 0
            while RANcorrectinput == 0:
                RANselection = input("Selection:")
                if RANselection == "1" or RANselection == "2":
                    RANcorrectinput=1
                else:
                    print("come on now, pick either 1 or 2")
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
                    if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                        if RANselection == "1":
                            Vmat[i,j] = ran.uniform(lower,upper)
                        if RANselection == "2":
                            Vmat[i,j] = ran.randint(lower,upper)
        if INPselection == "2":
            BOU = input("Boundary value:")
            BOUcorrect = 0
            while BOUcorrect == 0:
                try:
                    BOU = float(BOU)
                except ValueError:
                    print("Pick a number")
                    BOU = input("Boundary value:")
                else:
                    BOUcorrect = 1
            BOU = float(BOU)
            for i in range(0,dim):
                for j in range(0,dim):
                    if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                        Vmat[i,j] = BOU
        if INPselection == "3":
            for i in range(0,dim):
                for j in range(0,dim):
                    if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                        Vmat[i,j] = i/((j+1)*dim )
        if INPselection == "4":
            for i in range(0,dim):
                for j in range(0,dim):
                    if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                        Vmat[i,j] = m.sqrt((i**2 + j**2)/(2*dim*dim))
    if METHselection == "2":
        L = input("Plate Length (percent of matrix length taken up by plate):")
        Lcorrect = 0                      
        while Lcorrect == 0:
            try:
                L = int(L)
            except ValueError:
                print("Pick a number between 0 and 100")
                L = input("Plate Length (percent of matrix length taken up by plate):")
            else:
                if L > 100 or L <0:
                    print("Pick a number between 0 and 100")
                    L = input("Plate Length (percent of matrix length taken up by plate):")
                else:
                    Lcorrect = 1
        L = int(L)
        for i in range(int(dim*(0.5 - L/200)),int(dim*(0.5 + L/200))):
            Len = int(dim*((0.5 - L/200) - (0.5 + L/200)))
            ChargeDensity[i, int(dim*(3/8))] = 1/Len
            ChargeDensity[i, int(dim*(5/8))] = -1/Len                
    if METHselection != "3":
        #print(ChargeDensity, Vmat)
        ANS = Relax(ChargeDensity, Vmat, h, dim)
        plt.imshow(ANS, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.show()
    if METHselection == "3":
        print("DEFINE CONDITIONS")
        dim = input("Length of Matrix Dimensions:")
        DIMcorrect = 0
        while DIMcorrect == 0:
            try:
                DIMinput = int(dim)
            except ValueError:
                print("Pick an integer")
                dim = input("Length of Matrix Dimensions:")
            else:
                DIMcorrect = 1
        dim = int(dim)
        h = input("Length of Site (m):")
        Hcorrect = 0
        while Hcorrect == 0:
            try:
                Hinput = float(h)
            except ValueError:
                print("Pick an integer")
                h = input("Length of Site (m):")
            else:
                Hcorrect = 1
        h = float(h)
        tim = input("Total Time of Progression (s):")
        TIMcorrect = 0
        while TIMcorrect == 0:
            try:
                TIMinput = float(tim)
            except ValueError:
                print("Input a number")
                tim = input("Total Time of Progression (s):")
            else:
                TIMcorrect = 1
        tim = float(tim)
        tstep = input("Length of Time Step (s):")
        TSTEPcorrect = 0
        while TSTEPcorrect == 0:
            try:
                TSTEPinput = float(tstep)
            except ValueError:
                print("Input a number")
                tstep = input("Length of Time Step(s):")
            else:
                TSTEPcorrect = 1
        k = float(tstep)
        r = (k*alpha)/(h*h)
        rair = (k*alphaair)/(h*h)
        Tsteps = int(tim/k)
        Rod = 20*np.ones(dim)
        RodPrint = np.zeros((dim, Tsteps))
        ABSPrint = np.zeros(Tsteps)
        ABSPrint[0] = ABS(Rod)
        RodPrint[:,0] = Rod
        oven = input("Number of Sites in Oven:")
        OVEcorrect = 0
        while OVEcorrect == 0:
            try:
                OVEinput = int(oven)
            except ValueError:
                print("Input an integer")
                oven = input("Number of Sites in Oven:")
            else:
                OVEcorrect = 1
        oven = int(oven)
        ice = input("Number of Sites in Ice:")
        ICEcorrect = 0
        while ICEcorrect == 0:
            try:
                ICEinput = int(oven)
            except ValueError:
                print("Input an integer")
                ice = input("Number of Sites in Ice:")
            else:
                ICEcorrect = 1
        ice = int(ice)
        RotMat = RotMat(alpha, k, h, dim)
        OutPut = TimeProgress(tim, Rod, RotMat, alphaair, alphaice, k, Tsteps, ice, oven)[0]
        Display(dim, tim, k, OutPut, 0)
        Display(dim, tim, k, OutPut, 1)
    print("Finished?")
    print("[1] - Yes, quit")
    print("[2] - No, from the top")            
    QUITcorrectinput = 0 
    while QUITcorrectinput == 0:
        QUITselection = input("Selection:")
        if QUITselection == "1" or QUITselection == "2":            
            QUITcorrectinput=1
        else:
            print("come on now, pick 1 or 2")
    if QUITselection == "1":
        finish = 1