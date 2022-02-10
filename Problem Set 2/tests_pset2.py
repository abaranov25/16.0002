################################################################################
# 16.0002 / 18.0002 Fall 2021
# Problem Set 2: tests_pset2.py
# Name: Allen Baranov
# Collaborators:

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import IVP_numpy
from mysolve import mysolve
from mysolve import mysolve_vectorized


    
################################################################################
## test_IVP
# Solution to the system of ODE, time evolution
################################################################################

class Home_heating(IVP_numpy.IVP_numpy):
    
    def evalf(self, u, t):

        kb = self.get_p('kb')
        kA = self.get_p('kA')
        
        f = kA@u + kb
        
        return f
    
    
class Home_heating2(IVP_numpy.IVP_numpy):
    
    def evalf(self, u, t):
        
        kA = self.get_p('kA')
        kb = self.get_p('kb')(t)
        
        f = kA@u + kb
        
        return f


################################################################################
## test0
# Equilibrium solution
################################################################################

def run_test0():
        
    k = 0.5              # Newton cooling constant
    u_air = -3.0         # outside air temperature (Celsius)
    u_ground = 2.0       # near-surface ground temperature (Celsius)

    u_boundary = np.array([u_air, u_air, u_ground])
    
    power_per_room_heater = 10.0
    heaters_on = np.array([1.0,0.0,0.0])    # In order: annex, main, and basement
    power_per_room_GSHP = 6.0            # Ground-source heat pump
    COP_GSHP = 3.0                          # Coefficient of performance
    GSHP_on = np.array([0.0,1.0,0.0])       # In order: annex, main, and basement
    
    b_heaters = (power_per_room_heater*heaters_on)/k
    b_GSHP = (COP_GSHP*power_per_room_GSHP*GSHP_on)/k
    b = u_boundary + b_heaters + b_GSHP
    
    total_power = np.sum(power_per_room_heater*heaters_on + power_per_room_GSHP*GSHP_on)
    
    A = np.array( [[-2.0, 1.0, 0.0], [1.0, -3.0, 1.0], [0.0, 1.0, -2.0]] )
        
    u = mysolve(A,-b)

    s = ('Equilibrium temperature in: annex [' + format(u[0],'.2e') + 
         '] main ['+ format(u[1],'.2e') +
         '] basement ['  + format(u[2],'.2e') + '] ' +
         '     Total power: ' + format(total_power,'.2e'))

    return s
    


################################################################################
## test1
# Equilibrium solution, comparison between solvers
################################################################################
    

def run_test1():
       
    
    A = np.array( [[-2.0, 1.0, 0.0], [1.0, -3.0, 1.0], [0.0, 1.0, -2.0]] )

    b = np.array([20.0, 18.0, 0.0])
        
    u = linalg.solve(A,-b)
    v = mysolve(A,-b)
    w = mysolve_vectorized(A,-b)
    
    err_solve_vs_linalg = max(abs(v-u))
    err_solve_vectorized_vs_linalg = max(abs(w-u))

    s = ('Discrepancy mysolve vs linalg = ' + format(err_solve_vs_linalg,'.2e') +
         '         Discrepancy mysolve_vectorized vs linalg = ' + format(err_solve_vectorized_vs_linalg,'.2e'))

    return s




################################################################################
## test2
# Initial-value problem, time-varying forcing
################################################################################

def run_test2():
    
    
    tI = 0.0
    tF = 10.0            # hours
    
    k = 0.5              # Newton cooling constant (1/h)
    u_outside = -3.0     # outside air temperature (Celsius)
    u_ground = 2.0       # near-surface ground temperature (Celsius)
    
    uI = np.array([u_outside, u_outside, u_ground])    
    u_boundary = np.array([u_outside, u_outside, u_ground])

    power_per_room_heater = 10.0
    heaters_on = np.array([1.0,0.0,0.0])    # In order: annex, main, and basement
    power_per_room_GSHP = 6.0               # Ground-source heat pump
    COP_GSHP = 3.0                          # Coefficient of performance
    GSHP_on = np.array([0.0,1.0,0.0])       # In order: annex, main, and basement
    
    b_heaters = (power_per_room_heater*heaters_on)/k
    b_GSHP = (COP_GSHP*power_per_room_GSHP*GSHP_on)/k
    
    b = u_boundary + b_heaters + b_GSHP
    
    total_power = np.sum(power_per_room_heater*heaters_on + power_per_room_GSHP*GSHP_on)
    
    A = np.array( [[-2.0, 1.0, 0.0], [1.0, -3.0, 1.0], [0.0, 1.0, -2.0]] )
        
    p = {}
    p['kb'] = k*b
    p['kA'] = k*A
    
    this_IVP = Home_heating(uI, tI, tF, p)
    
    dt = 2e-1
    t, u = IVP_numpy.solve(this_IVP, dt, IVP_numpy.step_FE)
    
    plt.close('all')
    fig = plt.figure()
    plt.plot(t, u[:,0],'g-',label='Annex')
    plt.plot(t, u[:,1],'r-',label='Main')
    plt.plot(t, u[:,2],'b-',label='Basement')
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (deg C)')
    plt.grid(True)
    plt.legend() 
    
    

    
    ##### YOUR CODE GOES HERE #####
    
    p2 = {}
    
    def kb(t):
        if (t > 4):
            b2 = u_boundary + b_heaters + b_GSHP
        else:
            b2 = u_boundary + b_heaters
        
        return k*b2
    
    p2['kb'] = kb
    p2['kA'] = k*A
    
    this_IVP2 = Home_heating2(uI, tI, tF, p2)
    
    dt = 2e-1
    t2, u2 = IVP_numpy.solve(this_IVP2, dt, IVP_numpy.step_FE)
    
    fig = plt.figure()
    plt.plot(t2, u2[:,0],'g-',label='Annex')
    plt.plot(t2, u2[:,1],'r-',label='Main')
    plt.plot(t2, u2[:,2],'b-',label='Basement')
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (deg C)')
    plt.grid(True)
    plt.legend() 
    
    
       
    

################################################################################
## Main script
################################################################################

if __name__ == "__main__":
    
    
    # You may change the following list to set the tests you would like to run     
    tests_to_run = [run_test0, run_test1, run_test2]
    
    for test in tests_to_run:
        print()
        print(test.__name__,'results:')
        print(test())
            
    
