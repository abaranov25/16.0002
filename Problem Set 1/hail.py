################################################################################
# CSE.0002 Fall 2021
# Problem Set 1: hail 
# Name: Allen Baranov
# Collaborators:

import math
import matplotlib.pyplot as plt
import IVP

################################################################################
## Hail class definition (a subclass of IVP)
################################################################################

class HailIVP(IVP.IVP):
    def evalf(self, u, t):
        """
        Calculates forcing for hail trajectory given the current state and time.
        IMPORTANT: some of the provided helper functions assume u[0] is the velocity (V)
        and u[1] = the altitude (z).  And thus f[0] = dV/dt and f[1] = dz/dt
    
        Args:
            u (float list): current state.
            t (float): current time
    
        Returns:
            f (float list): current forcing. 
    
        """
        
        g    = self.get_p('g')
        rhoa = self.get_p('rhoa')
        rhop = self.get_p('rhop')
        Dp   = self.get_p('Dp')
        CD   = self.get_p('CD')
        
        pi = math.pi
        V = u[0]
    
        Ap = pi/4 * Dp**2
        D = 0.5*rhoa*(V**2)*Ap*CD
        m = pi/6*rhop*(Dp**3)
        
        f1 = g - D/m
        f2 = -V
        
        return [f1, f2]


################################################################################
## Functions to solve hail trajectory IVP
################################################################################

def hail_Verror(hail_IVP, t, V):
    """
    Calculates the exact (analytic solution) hail velocity Vex and the absolute error e    

    Args:
        hail_IVP (IVP object): IVP for hail trajectory
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        
    Returns:
        e (float list): velocity error e[n] = (abs(V[n]-Vex[n]))
        Vex (float list): exact (analytic solution) velocity values, Vex[n]

    """
    
    pi   = math.pi
    Dp   = hail_IVP.get_p('Dp')
    rhop = hail_IVP.get_p('rhop')
    g    = hail_IVP.get_p('g')
    CD   = hail_IVP.get_p('CD')
    rhoa = hail_IVP.get_p('rhoa')
    
    Ap = pi/4 * Dp**2
    m = pi/6*rhop*(Dp**3)
    Vterm = math.sqrt((2*m*g)/(rhoa*Ap*CD))
    
    dt = t[1]-t[0]
    V = IVP.solve(hail_IVP, dt, IVP.step_FE)
    
    e = []
    Vex = []
    
    for i in range(len(t)):
        Vex.append(Vterm * math.tanh(g*t[i]/Vterm))
        E = abs(Vex[i] - V[1][i][0])
        e.append(E)
        

    return e, Vex


def hail_Veplot(t, V, Vex, e, method='numerical'):
    """
    Produces a single figure with subplots of:
        * the numerical and exact velocity (overlayed) vs. time
        * the error vs. time
    The figure will be labeled with the method name

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        Vex (float list): exact (analytic solution) velocity values, Vex[n]
        e (float list): velocity error e[n]
        method (string, optional): name of the numerical method used to generate results

    """
    
    fig, (ax1, ax2) = plt.subplots(2,sharex=True)
    ax1.plot(t,Vex,"r-",label="exact")
    ax1.plot(t,V,"bo", label=method)
    ax2.plot(t,e, "bo")
    
    plt.xlabel("t (s)")
    ax1.set_ylabel("v (m/s)")
    ax2.set_ylabel("Error (m/s)")
    
    ax1.legend(loc="lower right")

    
def hail_Vzplot(t, V, z, method='numerical'):
    """
    Produces a single figure with subplots of:
        * velocity vs. time
        * z vs. time
    The figure will be labeled with the method name

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        z (float list): z values from numerical method, z[n]
        method (string, optional): name of the numerical method used to generate results

    """
    
    fig, (ax1, ax2) = plt.subplots(2,sharex=True)
    ax1.plot(t, V, "bo", label=method)
    ax2.plot(t, z, "bo")
    
    plt.xlabel("t (s)")
    ax1.set_ylabel("V (m/s)")
    ax2.set_ylabel("z (m)")
    
    ax1.legend(loc="lower right")


def hail_run_case(hail_IVP, dt, mlist):
    """
    Solves the hail trajectory problem described in hail_IVP with a timestep dt for 
    the methods in mlist and plots the solution and error as functions of time. 

    Args:
        hail_IVP (IVP object): Describes IVP case to be simulated
        dt (float): timestep
        mlist (function list): list of numerical methods to run on case

    """

    for method in mlist:
        t, u = IVP.solve(hail_IVP, dt, method)
        V = []
        z = []
        for i in range(len(t)):
            V.append(u[i][0])
            z.append(u[i][1])  
        e, Vex = hail_Verror(hail_IVP, t, V)
        emax = max(e)
        print('Max error = ',format(emax,'.3e'),' for ',method.__name__)
        hail_Veplot(t, V, Vex, e, method.__name__)                  
        hail_Vzplot(t, V, z, method.__name__) 
        
        
def hail_run_conv_study(hail_IVP, dtv, mlist, nflist):
    """
    Solves the hail trajectory problem described in hail_IVP for all of the timesteps
    in dtv and using all of the methods in mlist.  Then plots the error as a function 
    of dt and the number of function evaluations.

    Args:
        hail_IVP (IVP object): Describes IVP case to be simulated
        dtv (float list): timesteps 
        mlist (function list): list of numerical methods to run on case
        nflist (float list): for each method in mlist, gives the number of function
                evaluations pwer timestep the method uses

    """
    
    fig_evsdt, axs_evsdt = plt.subplots()
    fig_evsnf, axs_evsnf = plt.subplots()
    
    for m in range(len(mlist)):
        emax = []
        nf = []
        method = mlist[m]
        for dt in dtv:
            t, u = IVP.solve(hail_IVP, dt, method)
            V = []
            for i in range(len(t)):
                V.append(u[i][0])
            e, Vex = hail_Verror(hail_IVP, t, V)
            emax.append(max(e))
            nf.append(nflist[m]*len(t))
         
        axs_evsdt.loglog(dtv,emax,label=method.__name__)
        axs_evsnf.loglog(nf,emax,label=method.__name__)
        
    axs_evsdt.grid(True, which="both", ls="-")    
    axs_evsdt.set_xlabel('$\Delta t$')
    axs_evsdt.set_ylabel('Max error')
    axs_evsdt.legend(loc='lower right')
    
    axs_evsnf.grid(True, which="both", ls="-")    
    axs_evsnf.set_xlabel('Number of f evaluations')
    axs_evsnf.set_ylabel('Max error')
    axs_evsnf.legend(loc='lower left')
    
    

if __name__ == "__main__":
    uI = [0, 5000]
    tI = 0
    tF = 20 # s
    
    p = {}
    p['rhoa'] = 1.0 # kg/m^3
    p['rhop'] = 700.0 # kg/m^3
    p['Dp'] = 0.02 # m
    p['CD'] = 0.5
    p['g'] = 9.81 # m/s^2

    hail_IVP = HailIVP(uI, tI, tF, p)

# Set-up which methods to run
    mlist = [IVP.step_FE, IVP.step_RK4] # change this to only run one of the methods
    nflist = [1, 4] # number of f evaluations used per timestep each method

# Simulate a single case to demonstrate basic behavior of methods with dt=1.0
    hail_run_case(hail_IVP, 1.0, mlist)    
   
# Run varying dt for all methods and plot emax vs dt and emax vs function calls
    dtv = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]
    hail_run_conv_study(hail_IVP, dtv, mlist, nflist)
    

    
    
