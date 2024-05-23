# -*- coding: utf-8 -*-
"""
// ####################################################################################
// PDE Solver
// ####################################################################################
//
// File name   : GaussSeidel.py
// Version     : 1.0
// Author      : Di Lorenzo
// Date        : 23.05.2024
// Notes       : 
//
// ####################################################################################
"""

import matplotlib.pyplot as plt
import numpy as np


# Make the grid
xmin = 0
xmax = 2
Nx = 80
x,hx = np.linspace(xmin,xmax,Nx,retstep = True)
hx2 = hx**2
ymin = 0
ymax = 2
Ny = 40
y,hy = np.linspace(ymin,ymax,Ny,retstep = True)
hy2 = hy**2
X,Y = np.meshgrid(x,y,indexing='ij')
# Initialize potential
V = 0.5*np.ones_like(X)

# Dati e costanti

NA = 2e16 *10e6            # m^-3
ND = 1e16 *10e6            # m^-3
q = 1.602e-19         # C
epsilon0 = 8.85e-12   #F/m

# Enforce boundary conditions
V[:,0] = 0
V[:,-1] = 0
V[0,:] = 0
V[-1,:] = 0

# Allow possibility of charge distribution
#rho = np.zeros_like(X)
#rho = np.e**(X**2 - Y**2)
X11 = X[0:40, 0:40]
X12 = X[40:80, 0:40]

# RHO function
#rho_p = -q*NA*np.ones_like(X)  # C/cm^3   # Vettore X dalla prima colonna fino alla colonna Ny-1
#rho_n =  q*ND*np.ones_like(X)  # C/cm^3  # Vettore X dalla colonna Ny fino all'ultima colonna
#rho = 10e6*(rho_p + rho_n) / epsilon0   # 

rho_p = -q*NA*np.ones_like(X11)  # C/cm^3   
rho_n =  q*ND*np.ones_like(X12)  # C/cm^3  
# Concatenazione lungo l'asse delle colonne (verticale)
rho =  np.concatenate((rho_p, rho_n), axis=0) / epsilon0   # 




# Iterate
denom = 2/hx2 + 2/hy2
fig = plt.figure(1)
for n in range(200):
    # make plots every few steps
    if n % 10 == 0:
        plt.clf()
        ax = plt.axes(projection='3d')
        # Creating color map
        my_cmap = plt.get_cmap('hot')
        # Creating plot
        surf = ax.plot_surface(X,Y,V, cmap = my_cmap, edgecolor ='none')
        fig.colorbar(surf, ax = ax,
                     pad=0.15,
                     shrink = 0.3, 
                     aspect = 2)
        
        ax.set_title('Phi potential')
        #ax.set_zlim(-1, 2)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.draw()
        plt.pause(0.1)
# Iterate the solution
    for j in range(1,Nx-1):
        for k in range(1,Ny-1):
            V[j,k] = ( (V[j+1,k] + V[j-1,k])/hx2
                      +(V[j,k+1] + V[j,k-1])/hy2
                      +rho[j,k]) / denom