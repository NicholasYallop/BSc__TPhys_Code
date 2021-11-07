import math as m 
import numpy as np 
import random as ran
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import time
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D


def V(x, y, h):
    return((Vmat[x-1,y] + Vmat[x+1,y] + Vmat[x,y+1] + Vmat[x,y-1] + h*h*ChargeDensity[x,y])/4)
def Vcheck(x, y, h):
    return((Vmat[x-1,y] + Vmat[x+1,y] + Vmat[x,y-1] + Vmat[x,y+1] - 4*Vmat[x,y])/(h*h) + ChargeDensity[x,y])

            
Steps = 30
dim = 50
Length = int(dim*(6/8)) - int(dim*(3/8))
VPrint = np.zeros((Steps,Length))
TimePrint = np.zeros((Steps,2))
Iteration = np.zeros((Steps,2))
GradTest = np.zeros((Length, Steps))
GradTest2 = np.zeros((Length, Steps))
bound = 0
Test = np.zeros(Steps)
#for bound in range(0,2):
for aa in range(1,Steps):
    a = (0.4)*(aa/Steps)
    start_time = time.time()
    h = 0.02
    Vmat = np.zeros((dim, dim))
    VmatNew = np.zeros((dim, dim))
    Vcheckmat = np.zeros((dim-2, dim-2))
    ChargeDensity = np.zeros((dim, dim))
    """
    if bound == 0:
        for i in range(0,dim):
            for j in range(0,dim):
                if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                    Vmat[i,j] = i/((j+1)*dim)
        VmatNew[i,j] = Vmat[i,j]
    if bound == 1:
        for i in range(0,dim):
            for j in range(0,dim):
                if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                    Vmat[i,j] = m.sqrt((i**2 + j**2)/(2*dim*dim))
        VmatNew[i,j] = Vmat[i,j]
    """
    
    for i in range(int(dim*(0.5-a)),int(dim*(0.5+a))):
        L = int(dim*(0.5+a)) - int(dim*(0.5-a))
        ChargeDensity[i, int(dim*(3/8))] = 1/L
        ChargeDensity[i, int(dim*(6/8))] = -1/L
    """
    for i in range(0,dim):
        for j in range(0,dim):
            if i == dim-1 or j == dim-1 or i == 0 or j == 0:
                Vmat[i,j] = i/((j+1)*dim)
                VmatNew[i,j] = Vmat[i,j]
    """       
    #plt.imshow(Vmat, cmap='hot', interpolation='nearest')
    #plt.colorbar()
    #plt.show()
              
    n = 0
    Break = 0
    while Break == 0:
        Break = 1
        n += 1
        for i in range(1,dim-1):
            for j in range(1,dim-1):
                VmatNew[i,j] = V(i, j, h)
                if Vmat[i,j] != 0:
                    if abs((Vmat[i,j] - VmatNew[i,j])/Vmat[i,j]) > 10**-7:
                        Break = 0
                else:
                    if abs(Vmat[i,j] - VmatNew[i,j]) > 10**-8:
                        Break = 0
        Vmat = VmatNew
        #Vmat = 0.7*VmatNew + 0.3*Vmat
        if Break == 1:
            for i in range(1,dim-1):
                for j in range(1,dim-1):
                    HOLD = Vcheck(i,j,h)
                    if HOLD > 0.0001:
                        Break = 0
                    else:
                        Vcheckmat[i-1,j-1] = abs(HOLD)
    VDiff = Vmat[int(dim*0.5),int(dim*(3/8))]-Vmat[int(dim*0.5),int(dim*(6/8))]
    #ErrorPrint[aa, bound] = Vcheckmat.sum()/((dim-2)*(dim-2))
    VPrint[aa,:] = Vmat[int(dim*0.5), int(dim*(3/8)):int(dim*(6/8))]/VDiff
    for i in range(0,Length-1):
        GradTest[i,aa] = (VPrint[aa, i+1] - VPrint[aa, i])#/VDiff
        if i!=0:
            GradTest2[i,aa] = GradTest[i+1,aa] - GradTest[i,aa]
    Test[aa] = np.dot(GradTest2[:,aa],GradTest2[:,aa])/((Length-2)*(Length-2))
    TimePrint[aa, bound] = time.time()-start_time
    Iteration[aa, bound] = n
    print(aa, Vmat[int(dim*0.5),int(dim*(3/8))]-Vmat[int(dim*0.5),int(dim*(6/8))])

plt.plot(Vmat[int(dim*0.5), int(dim*(3/8)):int(dim*(6/8))])
plt.ylabel('Potential [V]')
plt.xlabel('No# Nodes from Left Plate')
plt.show()    
plt.plot(Test[5:Steps])
plt.ylabel('Sum of (d^2/dx^2)V across seperation Vm^-2')
plt.xlabel('Distance Between Plates [*0.013m]')
plt.show()
"""
plt.plot(ErrorPrint[:,0])
plt.plot(ErrorPrint[:,1])
plt.ylabel('Average Error')
plt.xlabel('Site Length (mm)')
plt.xlim((1,100))
plt.show()
"""
"""
plt.plot(TimePrint[:,0])
plt.plot(TimePrint[:,1])
plt.ylabel('Computation Time')
plt.xlabel('Site Length (mm)')
plt.xlim((1,100))
plt.show()
plt.plot(Iteration[:,0])
plt.plot(Iteration[:,1])
plt.ylabel('Number of Iterations')
plt.xlabel('Site Length ([e-3]m)')
plt.xlim((1,100))
plt.show()
"""
"""
plt.imshow(Vmat, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
#plt.imshow(Vcheckmat, cmap='hot', interpolation='nearest')
#plt.colorbar()
#plt.show()
"""


"""
plt.plot(Vmat[int(dim*0.5), int(dim*(2/5)):int(dim*(3/5))])
plt.ylabel('some numbers')
plt.show()
"""
fig = plt.figure()
ax = fig.gca(projection='3d')
# Make data.
X = np.arange(0, Length, 1)
Y = np.arange(0, Steps-5, 1)
X, Y = np.meshgrid(X, Y)
Z = VPrint[5:Steps,:]
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# Customize the z axis.
#ax.set_zlim(-0.0075, 0.0075)
ax.set_ylabel("Length of plate [*0.004m]")
ax.set_xlabel("No# Nodes from Left Plate")
ax.set_zlabel("Potential (Volts))")
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()