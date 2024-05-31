# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:38:52 2024

@author: dilorenzo
"""

import matplotlib.pyplot as plt
import numpy as np

##=====================================
# timer to count execution time
import time
start_time = time.time()
##=====================================

# Make the grid
xmin = 0
xmax = 50       # mm
Nx = 100
x,hx = np.linspace(xmin,xmax,Nx,retstep = True)
hx2 = hx**2
ymin = 0
ymax = 50        # m
Ny = 100
y,hy = np.linspace(ymin,ymax,Ny,retstep = True)
hy2 = hy**2
# =========================================================================
# Create necessary meshgrid and unravel grid from x,y. Define Nx, Ny
X,Y = np.meshgrid(x,y,indexing='ij')   # 2D meshgrid

#==========================================================================
# time vector
alpha = 10  # Thermal diffusivity
sim_time = 1 # seconds
dt = min(hx**2, hy**2) / (4 * alpha)
#dt = min((x[1]-x[0])**2, (y[1]-y[0])**2) / (4 * alpha)

# Initialize potential
V = np.full((Nx,Ny), 20.0) # Plate initially as 20 degres °C
#V[int(Nx/2), int(Ny/2)] = 20.0  # heating point
#V[2, 2] = 40.0  # heating point
#==============================================================================
# Enforce boundary conditions
# Dirichlet boundary conditions at outerwalls 
# (boundary condition type is defined through boundary operators)
V[:, 0] = np.linspace(0, 100, Ny)      # Left border (all the values of x and y=ymin)
V[:,-1] = np.linspace(0, 100, Ny)     # Right border (all the values of x and y=ymax)
V[0, :] = np.linspace(0, 100, Nx)      # Bottom border (all the values of y and x=xmin)
V[-1,:] = np.linspace(0, 100, Nx)     # Top border (all the values of y and x=xmax)


# Visualizing with a plot
fig, axis = plt.subplots()#plt.figure(1)
pcm = axis.pcolormesh(V, cmap=plt.cm.jet, vmin=0, vmax=100)
plt.colorbar(pcm, ax=axis)

# Simulating
counter = 0
denom = 2/hx2 + 2/hy2
while counter < sim_time:
        
    # Calcolo dell'errore tra le iterazioni successive
    previous_V = V.copy()  # Copia del potenziale precedente
    
    # Iterate the solution
    for j in range(1,Nx-1):
        for k in range(1,Ny-1):
            dd_ux = (V[j+1,k] + V[j-1,k] -2*V[j,k])/hx2
            dd_uy = (V[j,k+1] + V[j,k-1] -2*V[j,k])/hy2
            
            V[j,k] = alpha * dt * (dd_ux + dd_uy) + V[j,k]
            
    counter += dt
    print("t: {:.3f} [s], Average temperature: {:.2f} Celcius".format(counter, np.average(previous_V)))
    # Updating the plot
    pcm.set_array(previous_V)
    plt.xlabel('x')
    plt.ylabel('y')
    axis.set_title("Distribution at t: {:.3f} [s].".format(counter))
    plt.pause(0.01)
    #error_max = np.max(np.abs(V , previous_V)) # Calcolo dell'errore massimo
    #error= np.abs(V[j,k])/ error_max

plt.show()
##=====================================
# print the execution time    
print("---Simulation time = %s seconds ---" % (time.time() - start_time))

