################################################################################
# 16.0002 / 18.0002 Fall 2021
# Problem Set 3: climate modeling

import numpy as np
import math
import matplotlib.pyplot as plt
import IVP_numpy
from myrootfinding import bisection

################################################################################
## Climate class definition (a subclass of IVP_numpy)
################################################################################

class ClimateIVP(IVP_numpy.IVP_numpy):
    def evalf(self, u, t):
        """
        Calculates right-hand side in the climate model for atmospheric 
        temperature, and concentrations of H2O and CO2, given the current
        state and time

        Args:
            u (float list): current state.
            t (float): current time
    
        Returns:
            f (float list): current right-hand side. 
    
        """
        
        heatcapC = self.get_p('heatcapC')
        S0 = self.get_p('S0')
        sigma = self.get_p('sigma')
        tauH = self.get_p('tauH')
        tauC = self.get_p('tauC')
        H0 = self.get_p('H0')
        C0 = self.get_p('C0')
        PH = self.get_p('PH')
        PC = self.get_p('PC')
        alpha = self.get_p('alpha')
        eps = self.get_p('eps')
    
        T, H, C = u
        
        incfT = (1-alpha)*S0/4
        outfT = (1-eps(T,H,C)/2)*sigma*T**4
        
        fT = (1/heatcapC) * (incfT - outfT)
        fH = (1/tauH) * (H0 - H) + PH
        fC = (1/tauC) * (C0 - C) + PC
        f = np.array([fT, fH, fC])
        
        return f


################################################################################
## Test 0 : equilibrium temperature from root finding
################################################################################


def run_test0():

    def test_f(T):
        return T**2 - 1.0
    
    test_x1 = bisection(test_f,-2.0,0.1,60) + 1.0
    test_x2 = bisection(test_f,0.1,2.0,60) - 1.0
    test_bisection = max(abs(test_x1),abs(test_x2))

    p = {}
    p['heatcapC'] = 2.08e-8 # J/(m^2*K)
    p['S0'] = 1370 # W/m^2
    p['sigma'] = 5.67e-8 # W/(m^2*K^4)
    p['tauH'] = 2.5e-2 # yr
    p['tauC'] = 4e1 # yr
    p['H0'] = 0.85e4 # ppm
    p['C0'] = 280 # ppm
    p['PH'] = 0 # ppm
    p['PC'] = 5.0 # ppm/yr
    p['alpha'] = 0.3 
    p['eps'] = lambda T, H, C : 0.78

    Ta = 273 # K
    Tb = 300 # K

    climate_IVP = ClimateIVP(0, 0, 0, p)

    def f(T):

        u = np.array([T, p['H0'], p['C0']])
        rhs = climate_IVP.evalf(u,0)
        return rhs[0]

    x = bisection(f,Ta,Tb,10)
    
    s = ('run_test0: Test value = ' + format(test_bisection,'.3e') +
         '      Root near T = ' + format(x,'.3e'))
    print(s)


################################################################################
## Test 1 : time evolution
################################################################################

def run_test1(method):

    p = {}
    p['heatcapC'] = 2.08e-8 # J/(m^2*K)
    p['S0'] = 1370 # W/m^2
    p['sigma'] = 5.67e-8 # W/(m^2*K^4)
    p['tauH'] = 2.5e-2 # yr
    p['tauC'] = 4e1 # yr
    p['H0'] = 0.85e4 # ppm
    p['C0'] = 280 # ppm
    p['PH'] = 0 # ppm
    p['PC'] = 5.0 # ppm/yr
    p['alpha'] = 0.3     
    
    #p['eps'] = lambda T, H, C : 0.78 # reasonable value if one doesn't want to
                                      # use the function below
    def full_eps(T,H,C,C0):
        Patm = 1013.25  # hPa
        eps = 1.24*(H*Patm*(1e-6)/T)**(1/7) + 0.019/(math.log(2))*math.log(C/C0)
        return eps
    p['eps'] = lambda T, H, C: full_eps(T,H,C,p['C0'])
       
    TI = 287 # K
    HI = 0.85e4 # ppm
    CI = 820 # ppm
    uI = np.array([TI, HI, CI])

    tI = 0
    # Here choose tF appropriately for each question of the pset
    tF = 1e-6# yr

    climate_IVP = ClimateIVP(uI, tI, tF, p)
    
    # Here choose dt appropriately for each question of the pset
    dt = 1e-2  # be mindful of the number of steps needed to reach tF

    t, u = IVP_numpy.solve(climate_IVP, dt, method)

    plt.close('all')
    fig, axs = plt.subplots(3,1,sharex=True)
    axs_T = axs[0]
    axs_H = axs[1]
    axs_C = axs[2]
    axs_T.plot(t, u[:,0],'g-')
    axs_T.set_ylabel('T (K)')   
    axs_H.plot(t, u[:,1],'r-')
    axs_H.set_ylabel('H2O (ppm)')     
    axs_C.plot(t, u[:,2],'b-')    
    axs_C.set_ylabel('CO2 (ppm)')     
    axs_C.set_xlabel('t (yr)')
    plt.show()
    
if __name__ == "__main__":


# Set up which method to run
    
    # Here choose method appropriately for each question of the pset
    method = IVP_numpy.step_FE
    #method = IVP_numpy.step_BE
    # method = IVP_numpy.step_RK

# Comment out which test not to run
    run_test0()
    run_test1(method)
   

    
