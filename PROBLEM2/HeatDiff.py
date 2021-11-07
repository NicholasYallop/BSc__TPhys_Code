import math as m 
import numpy as np 
import random as ran
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


dim = 100
n = 20
l = 20
h = 0.01
k = 0.5
alpha = 59/(7900*450)
alphaair = 0.2435/(1006*1.276)
alphaice = 2.22/(916.2*2.05)
heattot = 0
r = (k*alpha)/(h*h)
rair = (k*alphaair)/(h*h)
T = 8000
Tsteps = int(T/k)
Rod = 20*np.ones(dim)
"""
for i in range(0,dim):
    Rod[i] =100*(1- (2/dim)**2 *(dim/2-i)**2)
"""

#Rod[int(dim/3):int(2*dim/3)-1] = 20*np.ones(int(dim/3))
"""
for i in range(0,dim):
    Rod[i] = (abs(12.5 - i))**1.15
    heattot += 1
"""
RodPrint = np.zeros((dim, Tsteps))
ABSprint = np.zeros(Tsteps)
ABSprint[0] = np.linalg.norm(Rod)
RodPrint[:,0] = Rod[0:dim]
Break = 0
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
           
RotMatInv = np.linalg.inv(RotMat)
for t in range(1,Tsteps):
    d = (1/k)*Rod    
    #print(np.linalg.det(RotMat))
    #print(RotMatInv)
    Rod = np.dot(RotMatInv, d)
    #print(Rod)
    """
    for i in range(0,n):
        Rod[i] = Rod[i] + (1000-Rod[i])*k*alphaair
    """
    for i in range(1,l):
        Rod[dim-i] = Rod[dim-i] + (-Rod[dim-i])*k*alphaice
    
    #print(dtMat)
    ABS = 0
    for i in range(0,dim):
        ABS += Rod[i]
    ABSprint[t]=ABS
    #print(ABS)
    RodPrint[:,t] = Rod[0:dim]
    print(t/Tsteps, ABS)

fig = plt.figure()
ax = fig.gca(projection='3d')
# Make data.
X = np.arange(0, dim, 1)
Y = np.arange(0, T, k)
X, Y = np.meshgrid(X, Y)
Z = np.transpose(RodPrint)
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# Customize the z axis.
#ax.set_zlim(0, 100)
ax.set_ylabel("Time (s)")
ax.set_xlabel("Distance Along Rod (m)")
ax.set_zlabel("Temperature (Degrees C)")
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

fig = plt.figure()
ay = fig.gca(projection='3d')
# Make data.
X = np.arange(0, dim, 1)
Y = np.arange(0, T, k)
X, Y = np.meshgrid(Y, X)
Z2 = RodPrint
# Plot the surface.
surf = ay.plot_surface(X, Y, Z2, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# Customize the z axis.
#ay.set_zlim(0, 100)
ay.set_xlabel("Time (s)")
ay.set_ylabel("Distance Along Rod (m)")
ay.set_zlabel("Temperature (Degrees C)")
ay.zaxis.set_major_locator(LinearLocator(10))
ay.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
"""
X=np.linspace(0,T, Tsteps)
plt.plot(X, ABSprint)
plt.show
"""

"""
    RodDisp = np.zeros((dim, 3))
    for i in range(0,3):
        RodDisp[:,i] = Rod[:]
    plt.imshow(RodDisp, cmap='hot', interpolation='nearest')
    plt.ylabel('some numbers')
    plt.show()
    """
"""
def dx(inpmat, i):
    if i!= 0 :
        dxminus = (inpmat[i] - inpmat[i-1])/h
    if i!= dim-1 :
        dxplus = (inpmat[i+1] - inpmat[i])/h
    if i != 0 and i!= dim-1:
        dx = 0.5*(dxminus + dxplus)
    elif i != 0:
        dx = dxminus
    else:
        dx = dxplus
    return(dx)
SpacialDiff = np.zeros(dim)
for i in range(0,dim):
    SpacialDiff[i] = dx(Rod, i)
"""