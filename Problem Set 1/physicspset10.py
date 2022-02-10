#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 20:22:17 2021

@author: allenbaranov
"""


import numpy as np
import matplotlib.pyplot as plt

"""
r1I,r2I,r3I denote the initial positions of the masses
v1I,v2I,v3I denote the initial velocities of the masses
G.M are physical constants provided by the problem
dt is the time step provided by the problem, in seconds
"""

r1I = np.array([-7.48e11, 1.496e11])
r2I = np.array([-7.48e11, -1.496e11])
r3I = np.array([1.496e12, 0])

v1I = np.array([-29788.2, 2371.89])
v2I = np.array([29788.2, 2371.89])
v3I = np.array([0, -4743.7])

G = 6.674e-11
M = 1.989e30
dt = 315360

"""
evalr and evalv compute the differential equations for r and v respectively

evalr returns an array containing [dr1/dt, dr2/dt, dr3/dt]
evalv returns an array containing [dv1/dt, dv2/dt, dv3/dt]
"""

def evalr(v1,v2,v3):
    return np.array([v1,v2,v3])

def evalv(r1,r2,r3):
    dv1 = G*M*((r3-r1)/np.linalg.norm(r3-r1)**3 + (r2-r1)/np.linalg.norm(r2-r1)**3)
    dv2 = G*M*((r3-r2)/np.linalg.norm(r3-r2)**3 + (r1-r2)/np.linalg.norm(r1-r2)**3)
    dv3 = G*M*((r1-r3)/np.linalg.norm(r1-r3)**3 + (r2-r3)/np.linalg.norm(r2-r3)**3)
  
    return np.array([dv1,dv2,dv3])

"""
stepr_FE and stepv_FE compute one step forward using the Forward-Euler method

stepr_FE returns an array containing the new position vectors
stepv_FE returns an array containing the new velocity vectors

"""

def stepr_FE(r1,r2,r3,v1,v2,v3):
    dr = evalr(v1,v2,v3)
    return np.array([r1,r2,r3]) + dt*dr

def stepv_FE(r1,r2,r3,v1,v2,v3):
    dv = evalv(r1,r2,r3)
    return np.array([v1,v2,v3]) + dt*dv

"""
This for loop finds the positions and velocities of the three masses and saves 
the position history for plotting

"""
N = 10000 # So that N * dt = 100 years
traj1 = np.zeros((N,2))
traj2 = np.zeros((N,2))
traj3 = np.zeros((N,2))

r1, r2, r3 = r1I, r2I, r3I
v1, v2, v3 = v1I, v2I, v3I

for i in range(N):
    traj1[i] = r1
    traj2[i] = r2
    traj3[i] = r3
    
    r = stepr_FE(r1,r2,r3,v1,v2,v3)
    v = stepv_FE(r1,r2,r3,v1,v2,v3)
    
    r1, r2, r3 = r[0],r[1],r[2]
    v1, v2, v3 = v[0],v[1],v[2]
    

"""
The following code converts from meters to astronomical units (AU) of length
"""

traj1 = traj1 / 1.496e11
traj2 = traj2 / 1.496e11
traj3 = traj3 / 1.496e11


"""
The following code plots the findings, and is provided in the problem
"""
fig, ax = plt.subplots(1,8, figsize=(16,2.0), sharey=True, sharex=True)
for i in range(8):
    max_t = i*(len(traj1)//8)
    ax[i].plot(traj1[:max_t,0], traj1[:max_t,1], label='M1')
    ax[i].plot(traj2[:max_t,0], traj2[:max_t,1], label='M2')
    ax[i].plot(traj3[:max_t,0], traj3[:max_t,1], label='M3')
    # add final points
    ax[i].plot(traj1[max_t,0], traj1[max_t,1], marker='o', color='tab:blue')
    ax[i].plot(traj2[max_t,0], traj2[max_t,1], marker='o', color='tab:orange')
    ax[i].plot(traj3[max_t,0], traj3[max_t,1], marker='o', color='tab:green')
    # add initial points
    ax[i].plot(traj1[0,0], traj1[1,1], marker='x', color='tab:blue')
    ax[i].plot(traj2[0,0], traj2[1,1], marker='x', color='tab:orange')
    ax[i].plot(traj3[0,0], traj3[1,1], marker='x', color='tab:green')
    ax[i].set_xlim(-10.0, 10.0)
    ax[i].set_ylim(-10.0, 10.0)
    ax[i].set_aspect(1)
plt.show()

"""
The following code converts from meters to astronomical units (AU) of length
"""
    