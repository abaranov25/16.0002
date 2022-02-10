# CSE.0002 Fall 2021
# Problem Set 4: *********
# Name: Allen Baranov
# Collaborators: 

import numpy as np
import matplotlib.pyplot as plt
import gd

def J0(xv):
    """
    Calculates J0 objective function and its gradient (see pset desription)

    Args:
        xv (1D numpy array): for this problem, only one value (xv[0]) will be present.

    Returns:
        J0 (float): objective function at x=xv[0].
        J0prime (1D numpy array): gradient of J0 at x=xv[0]

    COMMENT: YOU DO NOT NEED TO MODIFY THIS FUNCTION
    """
    x = xv[0]
    J0 = (x-0.5)**2
    J0prime = 2*(x-0.5)
    
    return J0, J0prime

def J1(xv, pdict):
    """
    Calculates J1 objective function and its gradient (see pset desription)

    Args:
        xv (1D numpy array): x=xv[0], y=xv[1]
        pdict (dictionary): contains key='user' with the associated value refering to the
            user location such that:
                pdict['user'][0] is the x-location of the user
                pdict['user'][1] is the y-location of the user

    Returns:
        J1 (float): objective function at xv
        J1prime (1D numpy array): gradient of J1 at xv
        
    """
    
    # Initializing variables
    xu = pdict['user'][0]
    yu = pdict['user'][1]
    
    x = xv[0]
    y = xv[1]
    
    dx = (x - xu)
    dy = (y - yu)
    
    # Finding J1
    J1 = -1/(1 + dx**2 + dy**2)
    
    # Finding gradients (use chain rule)
    J1primeX = J1**2 * 2 * (x - xu) 
    J1primeY = J1**2 * 2 * (y - yu) 
    
    J1prime = np.array([J1primeX, J1primeY])
    
    return J1, J1prime

def run_test0():
    # COMMENT: YOU DO NOT NEED TO MODIFY THIS FUNCTION
    
    print('Running test0')
    xstart = np.array([-0.5])
    alpha = 0.1
    nStop = 40
    xhist, Jhist = gd.gradient_descent(J0, xstart, alpha, nStop, verbose=True)
    print('In run_test0: xopt = {:.2e}, Jmin = {:.2e}'.format(xhist[-1][0],Jhist[-1]))

def run_test1():
    print()
    print('Running test1')
    pdict = {}
    pdict['user'] = np.array([0.5,-0.25])
    xstart = np.array([0.,0.])
    alpha = 0.1
    nStop = 40
    xhist, Jhist = gd.gradient_descent(J1, xstart, alpha, nStop, verbose=True, pdict=pdict)
    print('In run_test1: (xopt, yopt) = ({:.2e}, {:.2e}), Jmin = {:.2e}'.format(xhist[-1,0],xhist[-1,1],Jhist[-1]))
 
    # Set-up linearly spaced points in x and y for evaluating objective function
    Nx = 101
    Ny = 101
    bx = np.linspace(-1.,1.,Nx)
    by = np.linspace(-1.,1.,Ny)
    
    # Loop over all (x,y) points and evaluation objective function (call J1)
    
    J = np.zeros((Nx,Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            xv = [bx[i],by[j]]
            J[i,j] = J1(xv,pdict)[0]
    
    # Plot contours.  Include contour labels.  
    # Also set the axis to 'square' so the x and y axis will have the same scale (use axs.axis('square'))
    # Don't forget to label the x and y axes!  And turn on the grid
    
    # Apply a transpose to convert from row-column notation to x-y notation
    
    fig,axs = plt.subplots()
    c = axs.contour(bx,by,J.transpose())
    
    axs.clabel(c)
    axs.grid(True)
    axs.axis('square')
    
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks(np.arange(-1.,1.1,0.25))
    plt.yticks(np.arange(-1.,1.1,0.25))

    # Plot xhist markers on top of contours (use pyplot scatter function, i.e. axs.scatter(...)).  
    # Use green circles for the markers.    
    
    axs.scatter(xhist[:,0], xhist[:,1], color = "Green")

    # Plot final xhist marker on top of contours but use a magenta circle for the marker 
    
    axs.scatter(xhist[-1,0], xhist[-1,1], color = "Magenta")
        
    plt.show()
    
#run_test0()
run_test1() # Uncomment this when you are implementing J1
