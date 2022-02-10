# CSE.0002 Fall 2021
# Problem Set 4: *********
# Name: Allen Baranov
# Collaborators: 

import gd
import numpy as np
import matplotlib.pyplot as plt

def calc_objfun(basesv, pdict):
    """
        Calculates the objective function using the Lp-norm version of min and max to approximate 
        the min max power (min over all users, max over all bases).  See pset description for more
        details. 

    Args:
        basesv (1D numpy array): for base i, the (x,y) location = (basesv[2*i], basesv[2*i+1])
        pdict (dictionary) with the following key, val pairs:
            pdict['p'] = value of p used in Lp-norm min/max
            pdict['r0'] = r0 in power model (see pset description)
            pdict['users'] = 2D numpy array of user locations such that 
                user[i,0] is the x location of user i
                user[i,1] is the y location of user i

    Returns:
        Pp: objective function value
        Pp_prime: gradient of Pp with respect to basesv (same shape as basesv)
        
    COMMENT: YOU DO NOT NEED TO MODIFY THIS FUNCTION
    """
    p = pdict['p']
    r0 = pdict['r0']
    r02 = r0**2
    users = pdict['users']
    Nusers = len(users)
    Nbases = int(len(basesv)/2)
    bases = basesv.reshape((Nbases,2))
    Ppt = 0.
    Ppt_prime = np.zeros(bases.shape)
    Up_prime = np.zeros(bases.shape)
    for nu in range(Nusers):
        xu = users[nu,0]
        yu = users[nu,1]
        Up = 0.
        Up_prime[:,:] = 0.
        for nb in range(Nbases):
            xb = bases[nb,0]
            yb = bases[nb,1]
            d2ub = (xb-xu)**2 + (yb-yu)**2
            d2ub_xb = 2*(xb-xu)
            d2ub_yb = 2*(yb-yu)
            su = (r02 + d2ub)**(-p)
            su_d2ub = -p*(r02 + d2ub)**(-p-1)
            Up += su
            Up_prime[nb,0] +=  su_d2ub*d2ub_xb
            Up_prime[nb,1] +=  su_d2ub*d2ub_yb
            
        Ppt += (Up/Nbases)**(-1)
        Ppt_prime -= (Up/Nbases)**(-2.0)*Up_prime/Nbases
    Pp = -(Ppt/Nusers)**(-1./p)
    Pp_prime = (1.0/p)*(Ppt/Nusers)**(-1./p-1)*Ppt_prime/Nusers

    return Pp, Pp_prime.reshape(2*Nbases)


def plot_objfun_onebase(axs, pdict, Nx=101, Ny=101):
    """
    Plots objfun contours as a function of the (x,y) location of a single base.
    The (x,y) plotting region is chosen to be an additional 10% beyond the bounding box 
    of user locations.  In other words:
        xmin = minimum x location of all users
        xmax = maximum x location of all users
        ymin = minimum y location of all users
        ymax = maximum y location of all users  
        
        Lx = xmax - xmin
        Ly = ymax - ymin
        
        Then, the plotting region is for all points: 
            x0 < x < x1 and y0 < y < y1
        where:
            x0 = xmin - 0.1 Lx
            x1 = xmax + 0.1 Lx
            y0 = ymin - 0.1 Ly
            y1 = ymax + 0.1 Ly

    Overlayed on the plot are the user locations (this is done by calling plot_users)
    
    In addition, add the title to the plot giving the minimum value and the (x,y) location 
    of the minimum (use two digits beyond the decimal in the title)
    
    Args:
        axs: axs on which plotting will be done
        pdict (dictionary): dictionary as described in calc_objfun docstring.
        Nx (int, optional): Number of points in x for plotting (between x0 and x1) Defaults to 101.
        Ny (int, optional): Number of points in y for plotting (between y0 and y1) Defaults to 101.

    Returns:
        None.

    """
    
    users = pdict['users']
    
    # Find xmin, xmax, ymin, ymax
    
    x = users[:,0]
    y = users[:,1]
    
    xmin = x[np.argmin(x)]
    xmax = x[np.argmax(x)]
    ymin = y[np.argmin(y)]
    ymax = y[np.argmax(y)]    
    
    # Find Lx and Ly
    
    Lx = xmax - xmin
    Ly = ymax - ymin
    
    # Find x0, x1, y0, y1
    
    x0 = xmin - 0.1 * Lx
    x1 = xmax + 0.1 * Lx
    y0 = ymin - 0.1 * Ly
    y1 = ymax + 0.1 * Ly
    
    # Set-up linearly spaced points in x and y for evaluating objective function
    
    bx = np.linspace(x0,x1,Nx)
    by = np.linspace(y0,y1,Ny)
    
    # Loop over all (x,y) points and evaluation objective function (call calc_objfun)
    
    J = np.zeros((Nx,Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            xv = np.array([bx[i],by[j]])
            J[i,j] = calc_objfun(xv, pdict)[0]
            
    
    # Plot user location
    plot_users(axs, users)
    
    # Plot contours (on same axs).  Include contour labels.
    
    c = axs.contour(bx,by,J.transpose())
    axs.clabel(c)

    # Find minimum objective function and (x,y) location
    # TODO: Enter code here to do this.  Note: There are many ways to do this... 
    # but consider using np.argmin and np.unravel_index for a very numpy approach.
    
    minf = J[0,0]
    minfx = 0
    minfy = 0
    
    for i in range(Nx):
        for j in range(Ny):
            xv = np.array([bx[i],by[j]])
            current = calc_objfun(xv, pdict)[0]
            if current < minf: 
                minfx = bx[i]
                minfy = by[j]
                minf = current
    
    # Add a red asterick (*) at minimum
    
    axs.scatter(minfx, minfy, color = "Red", marker = "*")
    
    # Add title giving objective function minimum and (x,y) location in scientific
    # notation with 2 digits beyond the decimal.  
    
    title = "min f = {:.2e} at ({:.2e}, {:.2e})".format(minf,minfx,minfy)
    axs.set_title(title)


def plot_users(axs, users):
    """
    Plot the location of the users

    Args:
        axs: axs on which plotting will be done
        users (numpy 2D array): location of users (see calc_objfun docstring for more)
        
    COMMENT: YOU DO NOT NEED TO MODIFY THIS FUNCTION
    """
    axs.scatter(users[:,0],users[:,1],marker='s',color='blue')
    axs.set_xlabel('x')
    axs.set_ylabel('y')
    axs.axis('square')
    axs.grid(True)


def plot_optimize(axs, bases0, alpha, nStop, verbose=False, pdict={}):
    """
    Optimizes the base locations (calling the gradient descent algorithm) and then plots the
    optimization history of the base locations overlayed with the user locations.

    Args:
        axs: axs on which plotting will be done
        bases0 (2D numpy array): initial location of bases such that
            bases0[i,0] is the x location of base i
            bases0[i,1] is the y location of base i
            
        alpha (float): gradient descent step size
        nStop (int): number of iteration of gradient descent to run
        
        verbose (bollean, optional): verbose option for gd.radient_descent. Defaults to False.
        
        pdict (dictionary): dictionary as described in calc_objfun docstring.

    COMMENT: YOU DO NOT NEED TO MODIFY THIS FUNCTION
    """
    
    Nbases = len(bases0)
    bases0v = bases0.reshape(2*Nbases)
    bhist, Jhist = gd.gradient_descent(calc_objfun, bases0v, alpha, nStop, verbose, pdict)
    Jmin = Jhist[-1]
    bmin = bhist[-1,:]
    print('   Jmin   = {:.2e}'.format(Jmin))
    plot_users(axs, pdict['users'])
    for i in range(Nbases):        
        axs.scatter(bhist[ :,2*i], bhist[ :,2*i+1], color='green',  marker='o')
        axs.scatter(bmin[2*i], bmin[2*i+1], color='magenta', marker='o')  
        print('   Base {:d} = ({: 1.2e}, {: 1.2e})'.format(i,bmin[2*i],bmin[2*i+1]))


def optimize_bases123(users):
    """
    Performs optimization for 1, 2, and 3 base cases for the user locations in users.  For the 
    case of a single base, also call plot_objfun_onebase to plot the objective function.

    Args:
        users (2D numpy array): user location data (see docstring for calc_objfun)

    COMMENT: You do not need to make any change to this function except for adjustments to 
    alpha and Nstop.  No other changes are needed.
    """
    
    pdict = {}
    pdict['users'] = users
    pdict['r0'] = 1.0
    pdict['p'] = 7
    fig1, axs1 = plt.subplots()
    plot_objfun_onebase(axs1, pdict)

    # These values of alpha and nStop should be changed to improve the 
    # performance of gradient descents.  See pset description.
    alpha = 0.5
    nStop = 45
    
    print('Optimization with 1 base:')
    bases01 = np.array([[0.0,0.0]])
    plot_optimize(axs1, bases01, alpha, nStop, verbose=False, pdict=pdict)
    
    print()
    print('Optimization with 2 bases:')
    bases02 = np.array([[-0.5,-0.1],[0.1,0.5]])
    fig2, axs2 = plt.subplots()
    plot_optimize(axs2, bases02, alpha, nStop, verbose=False, pdict=pdict)

    print()
    print('Optimization with 3 bases:')
    bases03 = np.array([[-0.5,-0.1],[0.,0.],[0.1,0.5]])
    fig3, axs3 = plt.subplots()
    plot_optimize(axs3, bases03, alpha, nStop, verbose=False, pdict=pdict)  

# This is the main code below

# User locations for T-intersection study
users_T = np.array([[-1.,-1.],[-1,-0.5],[-1.,0.],[-1.0,0.5],[-1.0,1.0],[-0.5,0.],[0.,0.],[0.5,0.],[1.,0.]])

# User locations for undergrad residence study
users_UG = np.array([[-0.25,-0.35],[1.1,-0.12],[-0.5,-0.4],[-1.3,-0.4],[-1.0,0.35],[0.95,-0.12],[-1.8,0.2],[-0.65,1.1],[-1.65,-0.4],[-2.2,-0.4],[-1.86,-0.4],[-0.8,-0.4]])

optimize_bases123(users_UG) # Change this to users_UG or users_T depending on which case you want to run


