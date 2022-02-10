################################################################################
# CSE.0002 Spring 2021
# Problem Set 5: IVP_numpy


"""
A class representing an initial value problem (IVP), specifically of the form:
    
    du/dt = f(u,t,p)   for u(tI) = uI

In which the problem would be solved from t=tI to tF
And p are a set of parameters (fixed for all time)

Below the class definition are functions which implement numerical methods to approximate the solution to an IVP as described by the IVP class.

"""
import numpy as np
from scipy import optimize

class IVP_numpy():
    def __init__(self, uI, tI, tF, p = {}):
        """

        Args:
            uI (numpy ndarray): initial state.
            vI (numpy ndarray): initial time derivative of the state.
            tI (float): initial time.
            tF (float): final time.
            p (dictionary): set of fixed parameters. Optional.
        """
        self.uI = uI
        self.tI = tI
        self.tF = tF
        self.p  = p
        
    def evalf(self,u,t):
        """    

        Args:
            u (numpy ndarray): current state.
            t (float): current time.

        Returns:
            numpy ndarray: f(u,t).

        """
        raise NotImplementedError("evalf is not implemented for this object") 
    
    def get_tI(self):
        """        

        Returns:
            float: initial time.

        """
        return self.tI

    def get_tF(self):
        """

        Returns:
            float: final time.

        """
        return self.tF
    
    def get_uI(self):
        """        

        Returns:
            numpy ndarray: initial state.

        """
        return self.uI
    
    def get_p(self,name):
        """
      
        Arg:
            name (key): a key which should be in the object's parameter 
            dictionary 

        Returns:
            value of parameter key given by name

        """

        return self.p[name]
 

################################################################################
## Functions to numerically integrate an IVP
################################################################################

def step_FE(thisIVP,dt,un,tn):
    """
    Takes a single timestep of the Forward Euler method (FE) to 
    (approximately) integrate the state from u(tn) to u(tn+dt)
    
    Args:
        thisIVP (IVP_numpy object): object describing IVP being solved.
        dt (float): time increment.
        un (numpy ndarray): current state, i.e. u(tn).
        tn (float): current time.

    Returns:
        numpy ndarray: next state, i.e. u(tn+dt).

    """
    
    fn = thisIVP.evalf(un,tn)
    
    unplus1 = un + dt*fn
        
    return unplus1


def step_BE(thisIVP,dt,un,tn):
    """
    Takes a single timestep of the Backward Euler method (FE) to 
    (approximately) integrate the state from u(tn) to u(tn+dt)
    
    Args:
        thisIVP (IVP_numpy object): object describing IVP being solved.
        dt (float): time increment.
        un (numpy ndarray): current state, i.e. u(tn).
        tn (float): current time.

    Returns:
        numpy ndarray: next state, i.e. u(tn+dt).

    """
    def evalr_BE(u):
        return (u-un)/dt - thisIVP.evalf(u,tn+dt)
        
    un1 = optimize.root(evalr_BE, un)
    
    return un1.x

def step_RK4(thisIVP,dt,un,tn):
    """
    Takes a single timestep of the 4th order Runge-Kutta method (RK4) to
    (approximately) integrate the state from u(tn) to u(tn+dt)
    
    Args:
        thisIVP (IVP_numpy object): object describing IVP being solved.
        dt (float): time increment.
        un (numpy ndarray): current state, i.e. u(tn).
        tn (float): current time.

    Returns:
        numpy ndarray: next state, i.e. u(tn+dt).

    """
    
    k1 = thisIVP.evalf(un, tn)
    uns = un + 0.5*dt*k1
    k2 = thisIVP.evalf(uns, tn + 0.5*dt)  
    uns = un + 0.5*dt*k2
    k3 = thisIVP.evalf(uns, tn + 0.5*dt)
    uns = un + dt*k3
    k4 = thisIVP.evalf(uns, tn +     dt)
    uns = un + dt*(k1+2*k2+2*k3+k4)/6
            
    return uns



def solve(thisIVP, dt, method):
    """
    Solves an IVP using a timestep dt and method. Integrate from t=tI until u(tn) is 
    determined for which tn >= tF.

    Args:
        thisIVP (IVP_numpy object): object describing IVP being solved.
        dt (float): time increment.
        method (function): numerical integration method to use.

    Returns:
        t (numpy ndarray): time values at which u(t) is approximated. The nth item in 
            the list is the time of the nth step, tn = t[n].
        u (numpy ndarray): The values of the states at each step, so that 
            u[n,:] contains the state at time tn

    """

    # Sets initial conditions
    tI = thisIVP.get_tI()
    tF = thisIVP.get_tF()
    uI = thisIVP.get_uI()
    
    n = len(uI)
    t = np.arange(tI,tF+dt,dt)
    nt = len(t)
    
    u = np.zeros((nt, n))
    
    u[0,:] = uI
    
    # Loop from t=tI to t>=tF        
    for i in range(1,nt):
        u[i,:] = method(thisIVP,dt,u[i-1,:],t[i-1])
                
    return t, u



       
