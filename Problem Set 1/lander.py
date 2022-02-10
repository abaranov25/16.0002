################################################################################
# CSE.0002 Fall 2021
# Problem Set 1: lander
# Name: Allen Baranov
# Collaborators:

import math
import matplotlib.pyplot as plt
import IVP


################################################################################
## Lander class definition (a subclass of IVP)
################################################################################

class LanderIVP(IVP.IVP):

    def evalf(self, u, t):
        """
        Calculates forcing for a Martian lander given the current state and time.
    
        Args:
            u (float list): current state.
            t (float): current time
    
        Returns:
            f (float list): current forcing. 
    
        """
        
        m_l  = self.get_p('m_l')
    
        A_l  = self.get_p('A_l')
        CD_l = self.get_p('CD_l')
    
        A_p  = self.get_p('A_p')
        CD_p = self.get_p('CD_p')
        
        # You must use the lander.atmosphere method to get the atmospheric properties 
        # (the method has already been implemented for you immediately below)
    
        V = u[0]
        z = u[1]
        a = LanderIVP.atmosphere(self,z)
        rhoa = a.get("rhoa")
        g = a.get("g")
        
        Dl = 0.5*rhoa*(V**2)*A_l*CD_l
        Dp = 0.5*rhoa*(V**2)*A_p*CD_p
        
        pi = math.pi
        V_p = self.get_p('V_p')
        if (V > V_p):
            theta = self.get_p('theta_e')
            D = Dl
        else:
            theta = self.get_p('theta_p')
            D = Dl + Dp
        cosT = math.cos(pi/180 * theta)
        
        f = []
        f.append(g*cosT - D/m_l)
        f.append(-V*cosT)
        
        return f
        
        


    def atmosphere(self, z):
        """
        Returns properties of the Martian atmosphere at the altitude z (assumed in meters).
        The properties are returned using a dictionary.  Currently, the properties returned
        are:
            'Ta': temperature (K)
            'pa': pressure (kPa)
            'rhoa': density (kg/m^3)
            'g': gravity (m/s^2)
    
        Args:
            z (float): altitude in meters
    
        Returns:
            props (dictionary): property keys and values.  
    
        """
        
        props = {}
        if z>7000:
            Ta = 235 - 0.00108*(z-7000)
        else:
            Ta = 235 - 0.00222*(z-7000)
        props['Ta'] = Ta # K
            
        pa = 0.699*math.exp(-0.00009*z)
        props['pa'] = pa # kPa
         
        rhoa = pa/(0.1921*Ta) # kg/m^3
        props['rhoa'] = rhoa
        
        G = 6.674e-11 # Gravitational constant (m^3/ kg/s^2)
        Mm = 6.41693e23 # Mass of Mars (kg)
        Rm = 3389.5*1000 # Radius of Mars (m)
        
        g = G*Mm/(Rm+z)**2 # m/s^2
        props['g'] = g
        
        return props

################################################################################
## Functions to solve Martian lander trajectory IVP
################################################################################

def lander_Vzaplot(t, V, z, a, method='numerical'):
    """
    Produces a single figure with subplots of:
        * velocity vs. time
        * z vs. time
        * acceleration vs. time
    The figure will be labeled with the method name

    Args:
        t (float list): time values t[n]
        V (float list): velocity values from numerical method, V[n]
        z (float list): z values from numerical method, z[n]
        a (float list): acceleration values from numerical method, a[n]
        method (string, optional): name of the numerical method used to generate results

    """   
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)
    ax1.plot(t, V, "b-", label=method)
    ax2.plot(t, z, "b-")
    ax3.plot(t, a, "b-")
    
    ax1.set_ylabel("v (m/s)")
    ax2.set_ylabel("z (km)")
    ax3.set_ylabel("a (m/s^2)")
    plt.xlabel("t (s)")
   

def lander_run_case(lander_IVP, dt, m_list):
    """
    Solves the Martian lander problem described in lander_IVP with a timestep dt for 
    the methods in m_list.  For each method in m_list, this function will:
        * plot the solution by calling lander_Vzaplot.
        * print the time and altitude at which the parachute deployed.  
        * print the altitude and velocity at the final time

    Args:
        lander_IVP (IVP object): Describes IVP case to be simulated
        dt (float): timestep
        m_list (function list): list of numerical methods to run on case

    """
    V = []
    z = []
    a = []
    
    messagePrinted = False
    V_p = lander_IVP.get_p('V_p')
    
    
    for method in m_list:
        
        print("Method: " + method.__name__)
        print("-"*60)
        
        t, u = IVP.solve(lander_IVP, dt, method)
        for i in range(len(u)):
            V.append(u[i][0])
            z.append(u[i][1] / 1000)
            
            if ((not messagePrinted) and (V[i] < V_p)):
                print("Parachute deployed at t, z: " + str("{:.2e}".format(t[i])) + " s, " + str("{:.2e}".format( z[i])) + " km")
                messagePrinted = True
            
            a.append(LanderIVP.evalf(lander_IVP, u[i], t[i])[0])

        lander_Vzaplot(t, V, z, a, method)
        
        print("Final z, V: " + str("{:.2e}".format(z[-1])) + " km, " + str("{:.2e}".format(V[-1])) + " m/s")


if __name__ == "__main__":
    uI = [5800.0, 125*1000.0] # Initial velocity (m/s) and altitude (m)
    tI = 0
    tF = 300.0 # s
    
    p = {}
    p['m_l']     = 3300.0 # kg
    p['V_p']     = 470.0 # parachute deploy velocity (m/s)
    p['A_l']     = 15.9 # m^2
    p['CD_l']    = 1.7
    p['A_p']     = 201.0 # m^2
    p['CD_p']    = 1.2
    p['theta_e'] = 83.0# flight path angle during entry, in degrees
    p['theta_p'] = 70.0 # flight path angle during parachute, in degrees   
    

    lander_IVP = LanderIVP(uI, tI, tF, p)

# Set-up which methods to run
    m_list = [IVP.step_RK4] # change this to only run one of the methods

# Simulate a single case to demonstrate basic behavior of methods with dt=0.1 s
    lander_run_case(lander_IVP, 0.1, m_list)

    plt.show()
    

    
    
