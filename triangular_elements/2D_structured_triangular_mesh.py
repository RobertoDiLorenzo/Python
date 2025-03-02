# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 17:19:53 2024

@author: rober
"""

# https://www.youtube.com/watch?v=Pp__LuDcwoM&list=PLrVAK3EL234dN-N4Bnk1U2kKPEFsBeKtc&index=7&ab_channel=ComProgExpert

import numpy as np
import matplotlib.pyplot as plt

w = 10
h = 10

Nx = 3
Ny = 3

#allocate memory for mesh node
p = np.zeros((2, Nx*Ny))
# connectivity list
cl = np.zeros((3, 2*(Nx-1)*(Ny-1)))

#print(p)

dx = w / (Nx -1)
dy = h / (Ny -1)

index = 0
for i in range(0,Ny):
    
    y = (i)*dy
    for j in range(0,Nx):
        
        x = (j)*dx
        
        
        p[0, index] = x  # MATLAB uses 1-based indexing, Python uses 0-based indexing
        p[1, index] = y
        index += 1  # Increment index
        
print(p)

# Plot the mesh nodes
plt.plot(p[0, :], p[1, :], 'o')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Mesh Nodes')
plt.grid(True)
plt.show()

    

index = 0
for i in range(1, Ny):
    for j in range(1, Nx):
        index1 = j
        index2 = j + 1
        index3 = index2 + i * Nx
        index4 = index1 + i * Nx
        
        index += 1
        cl[:, index] = [index1, index3, index4]
        index += 1
        
        cl[:, index] = [index1, index2, index4]
        #index += 1

print(cl)