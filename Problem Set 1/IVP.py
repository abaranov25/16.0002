################################################################################
# CSE.0002 Fall 2021
# Problem Set 1: IVP
# Name: Allen Baranov
# Collaborators:

"""
A class representing an initial value problem (IVP), specifically of the form:
    
    du/dt = f(u,t,p)   for u(tI) = uI

In which the problem would be solved from t=tI to tF
And p are a set of parameters (fixed for all time)

Below the class definition are functions which implement numerical methods to approximate the solution to an IVP as described by the IVP class.

"""

class IVP():
    def __init__(self, uI, tI, tF, p = {}):
        """

        Args:
            uI (float list): initial condition of state.
            tI (float): initial time.
            tF (float): final time.
            p (dictionary): set of fixed parameters. Optional.
        """
        self.uI = uI[:]
        self.tI = tI
        self.tF = tF
        self.p  = p
        
    def evalf(self,u,t):
        """    

        Args:
            u (float list): current solution.
            t (float): current time.

        Returns:
            float list: f(u,t).

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
            float list: initial state.

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
        thisIVP (IVP object): object describing IVP being solved.
        dt (float): time increment.
        un (float list): current state, i.e. u(tn).
        tn (float): current time.

    Returns:
        float list: next state, i.e. u(tn+dt).

    """
    
    fn = thisIVP.evalf(un,tn)
    un1 = []
    for i in range(len(un)):
        un1.append(un[i] + dt*fn[i])
     
    return un1


def step_RK4(thisIVP,dt,un,tn):
    """
    Takes a single timestep of the 4th order Runge-Kutta method (RK4) to
    (approximately) integrate the state from u(tn) to u(tn+dt)
    
    Args:
        thisIVP (IVP object): object describing IVP being solved.
        dt (float): time increment.
        un (float list): current state, i.e. u(tn).
        tn (float): current time.
        

    Returns:
        float list: next state, i.e. u(tn+dt).

    """
    
    un1 = []
    fn = thisIVP.evalf(un,tn)
    
    ua = []; ub = []; uc = []; ud = [];
    fa = []; fb = []; fc = []; fd = [];
    
    #Finds first slope
    fa = fn
    for i in range(len(un)):
        ua.append(un[i] + 0.5*dt*fa[i])
        
    #Finds second slope
    fb = thisIVP.evalf(ua,tn+0.5*dt)
    for i in range(len(un)):
        ub.append(un[i] + 0.5*dt*fb[i])
        
    #Finds third slope
    fc = thisIVP.evalf(ub,tn+0.5*dt)
    for i in range(len(un)):
        uc.append(un[i] + dt*fc[i])
    
    #Finds fourth slope
    fd = thisIVP.evalf(uc,tn+dt)
    
    #Weighed Sum
    for i in range(len(un)):
        s = 1/6 * (fa[i] + 2*fb[i] + 2*fc[i] + fd[i])
        un1.append(un[i]+dt*s)
    
    return un1


def solve(thisIVP, dt, method):
    """
    Solves an IVP using a timestep dt and method. Integrate from t=tI until u(tn) is 
    determined for which tn >= tF.

    Args:
        thisIVP (IVP object): object describing IVP being solved.
        dt (float): time increment.
        method (function): numerical integration method to use.

    Returns:
        t (float list): time values at which u(t) is approximated. The nth item in 
            the list is the time of the nth step, tn = t[n].
        u (list of float lists): The values of the states at each step.  The nth 
            item in the list is the values of the states at tn.  i.e. u(tn) = u[n]
            where u[n] is a float list.  So, if there are three equations being integrated, then 
            u[n][0], u[n][1], and u[n][2] are the values of the three states at time t=t[n]

    IMPORTANT: The first element in the returned t and u lists will be the initial values 
    of t and u.  Thus:
        * t[0] will be a float which is equal to thisIVP.get_tI()
        * u[0] will be a float list which is equal to thusIVP.get_uI()
    """

    # Sets initial condition
    tI = thisIVP.get_tI()
    t = [tI]

    uI = thisIVP.get_uI()
    u = [uI]
    
    un = thisIVP.get_uI()
    tn = tI
    tF = thisIVP.get_tF()


    # Loop from t=tI to t>=tF
    while (tn < tF):
        #Using the provided method
        un1 = method(thisIVP, dt, un, tn)
        
        #Append values to outputs u and t
        u.append(un1)
        t.append(tn+dt)
        
        #Update un and tn for the next iteration
        un = un1
        tn = tn+dt
        
    return t, u




       
