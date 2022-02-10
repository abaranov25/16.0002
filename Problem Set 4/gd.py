# CSE.0002 Fall 2021
# Problem Set 4: *********
# Name: Allen Baranov
# Collaborators: 

import numpy as np

def gradient_descent(Jfunc, xstart, alpha, nStop, verbose=False, pdict={}):
    """
    Implements gradient descent algorithm to reduce an unconstrained objective function

    Args:
        Jfunc (function reference): reference to function that calculates the objective function and gradient.
            Jfunc is called with one required argument, x, the state at which the objective function and gradient
            are to be evaluated.  In addition, Jfunc will be passed pdict (if pdict is provided).  So, the return
            and call from a Jfunc will be:
                J, Jprime = Jfunc(x, pdict)    # if pdict was provided
                J, Jprime = Jfunc(x) # if pdict was not provided
                
        xstart (1D numpy array): Initial guess for x
        alpha (float): steepest descent step size factor
        nStop (int): number of iterations of steepest descent to run
        verbose (boolean): If true, print information at every iteration of the gradient descent.
        pdict (dictionary): passed through to Jfunc if provided.

    Returns:
        xhist (2D numpy array): returns the iteration history of x with xhist[n,:] corresponding 
                                to the nth iteration 
        Jhist (1D numpy array): returns the iteration history of the objective function
    """
    
    # Set up storage for xhist and Jhist
    Nx = len(xstart)
    xhist = np.zeros((nStop, Nx))
    Jhist = np.zeros(nStop)

    # Calculate J and Jprime at xstart
    x = xstart
    if len(pdict):
        J, Jprime = Jfunc(x, pdict)
    else:
        J, Jprime = Jfunc(x)
    xhist[0,:] = x
    Jhist[0] = J
    if verbose:
        print('Iter={:d}, J = {:.2E}'.format(0,J))
        
    # TODO: Complete the code inside this loop which will calculate for remaining iterations 
    # of the gradient descent algorithm
    for n in range(1,nStop):
        
        # Implement gradient descent here
        # Don't forget to output the required data if verbose is True.  See pset for more info.

        # There are many ways to find where the maximum component changes...
        # but consider using np.argmin and related functions for a very numpy approach.
        
        x1 = x - alpha * Jprime
        if len(pdict):
            J1, Jprime = Jfunc(x1, pdict)
        else:
            J1, Jprime = Jfunc(x1)
            
        dJ = J1 - J
        state = np.argmax(x1-x)
        dx = x1[state] - x[state]
        
        # Example of desired output
        # Iter=1: J = 6.40E-01, dJ = -3.60E-01, max dx = 2.00E-01 for state i=0
        if (verbose):
            print("Iter={:d}: J = {:.2e}, dJ = {:.2e}, max dx = {:.2e} for state i={:d}".format(n,J1,dJ,dx,state))
        
        
        #Saving values into output arrays
        xhist[n,:] = x1
        Jhist[n] = J1
        
        #Updating values for next iteration
        x,J = x1,J1
        
        
    return xhist, Jhist

