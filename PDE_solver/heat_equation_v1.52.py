# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:05:27 2024
// ####################################################################################
// PDE Solver
// ####################################################################################
//
// File name   : Heat_equation.py
// Version     : 1.52
// Author      : Di Lorenzo
// Date        : 7.06.2024
// Notes       : 
//
// ####################################################################################
"""

import matplotlib.pyplot as plt
import numpy as np
import time

import sys

def progressbar(prefix="", size=60, out=sys.stdout): 
    def show(j, count, start):
        x = int(size*j/count)
        remaining = ((time.time() - start) / j) * (count - j)        
        mins, sec = divmod(remaining, 60) 
        time_str = f"{int(mins):02}:{sec:03.1f}"
        out.write(f"{prefix}[{'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}\r")
        out.flush()
    
    sim_time = 100  # Esempio di tempo massimo
    counter = 0
    start_time = time.time()
    while counter < sim_time:
        # Esegui le operazioni del loop qui
        time.sleep(0.1)  # Simulazione di un'operazione
        
        counter += 1
        show(counter, sim_time, start_time)

    out.write("\n")
    out.flush()  
    
#================================================'
def square_source(j,k,heating_point_x,heating_point_y,heating_size):
    if (j>(int(heating_point_x) -heating_size) and j<(int(heating_point_x)+heating_size)) and (k>(int(heating_point_y) -heating_size) and (k<(int(heating_point_y) +heating_size))):
        state = True
    else:
        state = False
    return state

def circle_source(j,k,heating_point_x,heating_point_y,heating_size):
    """
    (x-xc)**2+(y-yc)**2 < r**2:
    """
    if ((j-int(heating_point_x))**2+(k-int(heating_point_y))**2 <= heating_size**2):
        state = True
    else:
        state = False
    return state
    

#
##=====================================
# timer to count execution time
start_time = time.time()
##=====================================
# Make the grid
# length
Lx = 5 #m
Ly = 5 #m
density_grid = 100   # [ptr/m^2]
ptr_grid= int(Lx*Ly*density_grid)


# Traslation
xmin = 0
xmax = Lx - xmin  # mm
ymin = 0
ymax = Ly - ymin  # m


Nx = int(Ly*density_grid)
x, hx = np.linspace(xmin, xmax, Nx, retstep=True)
hx2 = hx**2
Ny = int(Lx*density_grid)
y, hy = np.linspace(ymin, ymax, Ny, retstep=True)
hy2 = hy**2
# =========================================================================
# Create necessary meshgrid and unravel grid from x,y. Define Nx, Ny
X, Y = np.meshgrid(x, y, indexing='ij')   # 2D meshgrid

#==========================================================================
# Thermal Parameter
density = 2330  # density [kg/m^3]
cp = 700  # Specific heat [J/(kg*°C)]
k0 = 150  # thermal conductivity [J/s*m*°C]
alpha = k0/(density * cp) # Thermal diffusivity [m^2/s]   0.8[cm^2/s]

DEBUG=0

#==========================================================================
# time vector
sim_time = 5  # seconds
dt = min(hx**2, hy**2) / (4 * alpha)
#==========================================================================
# Initialize potential
V = np.full((Nx, Ny), 20.0)  # Plate initially as 20 degrees °C

# Initialize potential with heating point at the center
heating_point_x = int(Nx / 2)  # Coordinata x del punto di riscaldamento
heating_point_y = int(Ny / 2)  # Coordinata y del punto di riscaldamento
#Initialize source
S = np.full((Nx, Ny), 20.0)
heating_size = 1*density_grid
#==========================================================================


# Enforce boundary conditions
top_border    = 0#np.linspace(0, 100, Ny)      # Top border (all the values of ymax between xmin and xmax)
bottom_border = 0#np.linspace(0, 100, Ny)      # Bottom border (all the values of ymin between xmin and xmax)
left_border   = 0# np.linspace(0, 100, Nx)          # Left border (the values of xmin between ymin and ymax)
right_border  = 0#np.linspace(0, 100, Nx)      # Right border (the values of xmax between ymin and ymax)
# Dirichlet boundary conditions at outerwalls
# (boundary condition type is defined through boundary operators)
V[:, 0] = left_border        
V[:,-1] = right_border    
V[0, :] = top_border
V[-1,:] = bottom_border     
 
if DEBUG==1:
    # Visualizing with a plot
    #fig, axis = plt.subplots(1, figsize=(12, 6))  # Creazione di una figura con due assi (subplot)
    fig, axis = plt.subplots(1)
    
    # Grafico lungo Ly
    pcm1 = axis.imshow(V, cmap=plt.cm.jet, vmin=0, vmax=100, extent=[x.min(), x.max(), y.min(), y.max()])
    plt.colorbar(pcm1, ax=axis)
    axis.set_xlabel('x')
    axis.set_ylabel('y')


# Simulating
n = 0
counter = 0
#denom = 2/hx2 + 2/hy2
while counter < sim_time:
    previous_V = V.copy()  # Copia del potenziale precedente
    # Iterate the solution
    for j in range(1, Nx-1):
        for k in range(1, Ny-1):
            # Aggiungere il contributo della sorgente termica con il ciclo if
            if circle_source(j,k,heating_point_x,heating_point_y,heating_size):#square_source(j,k,Nx,Ny,heating_size):
                #if DEBUG == 1:
                    #print("non faccio nulla")
                S[j,k] = 100  # Celsius
                V[j,k]=S[j,k]
            else:
                dd_ux = (V[j+1, k] + V[j-1, k] - 2*V[j, k]) / hx2
                dd_uy = (V[j, k+1] + V[j, k-1] - 2*V[j, k]) / hy2
                V[j, k] = alpha * dt * (dd_ux + dd_uy) + V[j, k] 
                
                    
    counter += dt
    print("t: {:.3f} [s], Average temperature: {:.2f} Celcius".format(counter, np.average(V)))
    if DEBUG == 1:
        if n % 2 == 0:
            # Updating the plot
            pcm1.set_array(V)
        
            axis.set_title("Distribution at t: {:.3f} [s].".format(counter))
            plt.draw()
            plt.pause(0.02)
    n += 1

if DEBUG==0:
    # Visualizing with a plot
    #fig, axis = plt.subplots(1, figsize=(12, 6))  # Creazione di una figura con due assi (subplot)
    fig, axis = plt.subplots(figsize=(12, 6))
    # Grafico lungo Ly
    pcm1 = axis.imshow(V, cmap=plt.cm.jet, vmin=0, vmax=100, extent=[x.min(), x.max(), y.min(), y.max()])
    plt.colorbar(pcm1, ax=axis)
    axis.set_xlabel('x')
    axis.set_ylabel('y')
    
plt.tight_layout()
plt.show()
            
##=====================================
# print the execution time    
print("---Simulation time = %s seconds ---" % (time.time() - start_time))
