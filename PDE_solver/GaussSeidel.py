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
Nx = 160
x,hx = np.linspace(xmin,xmax,Nx,retstep = True)
hx2 = hx**2
ymin = 0
ymax = 2
Ny = 80
y,hy = np.linspace(ymin,ymax,Ny,retstep = True)
hy2 = hy**2
# =========================================================================
# Create necessary meshgrid and unravel grid from x,y. Define Nx, Ny
X,Y = np.meshgrid(x,y,indexing='ij')   # 2D meshgrid
# Initialize potential
V = 0.5*np.ones_like(X)

# Dati e costanti

NA = 2e16 *10e6            # m^-3
ND = 1e16 *10e6            # m^-3
q = 1.602e-19         # C
epsilon0 = 8.85e-12   #F/m


#==============================================================================
# Enforce boundary conditions
# Dirichlet/Neumann boundary conditions at outerwalls 
# (boundary condition type is defined through boundary operators)
V[:, 0] = 0      # Left border (all the values of x and y=ymin)
V[:,-1] = 0     # Right border (all the values of x and y=ymax)
V[0, :] = 0      # Bottom border (all the values of y and x=xmin)
V[-1,:] = 0     # Top border (all the values of y and x=xmax)

# Allow possibility of charge distribution
#rho = np.zeros_like(X)
#rho = np.e**(X**2 - Y**2)
X11 = X[0:Nx//2, 0:Ny]   # Nx//2 Divisione intera
X12 = X[Nx//2:Nx, 0:Ny]

# RHO function
rho_p = -q*NA*np.ones_like(X11)  # C/cm^3   
rho_n =  q*ND*np.ones_like(X12)  # C/cm^3  
# Concatenazione lungo l'asse delle colonne (verticale)
rho = np.concatenate((rho_p, rho_n), axis=0) / epsilon0   # 


# Iterate
max_iterations=200
denom = 2/hx2 + 2/hy2
fig = plt.figure(1)

for n in range(max_iterations):
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
                     shrink = 0.5, 
                     aspect = 2)
        
        ax.set_title('Phi potential')
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