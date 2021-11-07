import cmath
import math
import numpy as np
import random as rand
dim=2
iterations = 100
attempts = 100
Input=np.zeros((dim,dim))
"define random input for tests:"
for i in range(0,dim):
    for j in range(0,dim):
        Input[i,j]=rand.randint(1,10)
print(np.linalg.eig(Input)[1])
"eigenfinder"
eigenvects = np.zeros((dim,dim))
x = np.ones((dim, iterations))
eigencounter = 0
breakcounter = 0
for q in range(0, attempts):
    if breakcounter == 1:
        break
    for k in range(0,iterations):    
        if k ==0:
            for i in range(0,dim):
                x[i,0] = np.random.uniform(-1,1)
            diff = 0
            for r in range(0,dim): 
                for p in range(0,iterations):
                    x[0:dim,0] = x[0:dim,0] - eigenvects[0:dim,r]*np.dot(x[0:dim,0],eigenvects[0:dim,r])
            #for r in range(0,dim):
             #   print(eigenvects[0:dim,r], np.dot(x[0:dim,0],eigenvects[0:dim,r]))
            x[0:dim,0] = x[0:dim,0]/math.sqrt(np.dot(x[0:dim,0],x[0:dim,0]))            
        else:     
            stepup = np.dot(Input, x[0:dim,k-1])
            x[0:dim,k] = stepup/math.sqrt(np.dot(stepup,stepup))
            for r in range(0,dim):
                x[0:dim,k] = x[0:dim,k] - eigenvects[0:dim,r]*np.dot(x[0:dim,k],eigenvects[0:dim,r])
            x[0:dim,k] = x[0:dim,k]/math.sqrt(np.dot(x[0:dim,k],x[0:dim,k]))
    #print(x[0:dim, iterations-1] - eigenvects[0:dim,0], x[0:dim, iterations-1] + eigenvects[0:dim,0])
    eigenunique = 0
    for s in range(0,dim):
        if abs(np.dot(eigenvects[0:dim,s],x[0:dim, iterations-1])) <10**-3 and eigenvects[0:dim,s].any != x[0:dim, iterations-1].any and eigenvects[0:dim,s].any != (-1*x[0:dim, iterations-1]).any:
            eigenunique += 1
    if eigenunique == dim:
        eigenvects[0:dim,eigencounter] = x[0:dim, iterations-1] 
        eigencounter += 1
        if eigencounter == dim:
            print("done")
            breakcounter = 1
            print(eigenvects)
