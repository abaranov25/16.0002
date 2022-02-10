################################################################################
# CSE.0002 Fall 2021
# Problem Set 5: viral spread SEIR modeling
# Name:
# Collaborators:

import numpy as np
import math
import matplotlib.pyplot as plt
import IVP_numpy

################################################################################
## SEIR class definition (a subclass of IVP_numpy)
################################################################################

class SEIR_IVP(IVP_numpy.IVP_numpy):
    
    def evalf(self, u, t):
        """
        Calculates right-hand side in the SEIR model

        Args:
            u (numpy ndarray): current state.
            t (float): current time
    
        Returns:
            f (numpy ndarray): current right-hand side. 
    
        """
    
        calR = self.get_p('calR0')    # R0
        v0 = self.get_p('v0')         # vaccine effectiveness [0,1]
        Tv = self.get_p('Tv')         # decay time constant of v (days)
        Tinf = self.get_p('Tinf')     # infection time (days)
        Tinc = self.get_p('Tinc')     # incubation time (days)
        N = self.get_p('N')           # total population size
        vdecay = self.get_p('vdecay') # whether vaccine effectiveness decays
    
        S1, S2, E1, E2, I1, I2, R1, R2 = u

        vcurrent = v0

        # Question B.3        
        # Here compute vcurrent, to be used below as the vaccine efficiency
        pass  # TO DO... Remove this line when you implement this function

            
        fS1 = - calR * (1. - vcurrent) * (I1 + I2) * S1/ Tinf / N
        fS2 = - calR * (I1 + I2) * S2/ Tinf / N
        fE1 = (calR * (1. - vcurrent) * (I1 + I2) * S1/ Tinf / N) - E1/Tinc
        fE2 = (calR * (I1+I2) * S2/ Tinf / N) - E2/Tinc
        fI1 = E1/Tinc - I1/Tinf
        fI2 = E2/Tinc - I2/Tinf
        fR1 = I1/Tinf
        fR2 = I2/Tinf
        
        f = np.array([fS1, fS2, fE1, fE2, fI1, fI2, fR1, fR2])
        
        return f

        
################################################################################
## Single run of the SEIR model
################################################################################


def single_run(calR0, v0, beta, vdecay=False):
    """
    Computes the time evolution of the SEIR model

    Args:
        calR0, v0: random parameters in the ODE model (floats)
        beta: random parameter describing post-processing of ODE trajectory (float)
        vdecay: whether vaccine effectiveness decays in time (True/False) 

    """
            
    p = {}
    
    p['calR0'] = calR0       # reproductive number of virus
    p['v0'] = v0          # initial vaccine effectiveness
    p['beta'] = beta      # hospitalization rate
    p['vdecay'] = vdecay  # whether vaccine effect decays or not
    
    p['Tv'] = 200          # decay time constant of vaccine effectiveness
    p['Tinf'] = 3.0     # days
    p['Tinc'] = 5.0     # days
        
    p['N'] = 1e6

    # the '1' population is vaccinated and the '2' population is not
    EI1 = 700
    EI2 = 300
    SI1 = 700000
    SI2 = 299000

    II1 = 0
    II2 = 0
    RI1 = 0
    RI2 = 0

    uI = np.array([SI1, SI2, EI1, EI2, II1, II2, RI1, RI2])

    tI = 0
    tF = 3e2    # days

    this_SEIR_IVP = SEIR_IVP(uI, tI, tF, p)
    
    dt = 1
    
    t, u = IVP_numpy.solve(this_SEIR_IVP, dt, IVP_numpy.step_FE)
    
    return t, u

################################################################################
## Single run of the SEIR model, plot
################################################################################

def plot_single_run():
    
    calR0 = 4.0
    v0 = 0.75
    beta = 0.05

    t, u = single_run(calR0, v0, beta, vdecay=False)
    
    plt.close('all')
    fig = plt.figure()
    plt.plot(t, u[:,0],'g-',label='S vax')
    plt.plot(t, u[:,1],'g--',label='S no vax')
    plt.xlabel('t (days)')
    plt.ylabel('count (persons)')
    plt.plot(t, u[:,2],'r-',label='E vax')
    plt.plot(t, u[:,3],'r--',label='E no vax')

    plt.plot(t, u[:,4],'b-',label='I vax')
    plt.plot(t, u[:,5],'b--',label='I no vax')

    plt.plot(t, u[:,6],'k-',label='R vax')
    plt.plot(t, u[:,7],'k--',label='R no vax')

    plt.grid(alpha=0.25)
    plt.legend(loc='center right')
    plt.savefig('SEIR_single_run.pdf')
    plt.show()
    

################################################################################
## Statistics of the SEIR model with random parameters
################################################################################

    
def stats_many_runs(nruns):
    """
    Computes estimators for
        - the probability that the max of H is above Htot
        - the expectation of R1 at t=tF

    Args:
        nruns (int): number of times the SEIR model is run with independent random parameters

    """
    import random
    rng = np.random.default_rng()
    
    Hmax = np.zeros(nruns)
    R1tF = np.zeros(nruns)

    Htot = 1.5e3
    
    zstar = 1.96    # for 95% confidence intervals
    
    for i in range(nruns):

        # Here draw random numbers for calR0, beta, v0
        pass  # TO DO... Remove this line when you implement this function

        t, u = single_run(calR0, v0, beta, vdecay=False)
                
        # Question B.1, here compute Hmax[i]
        pass  # TO DO... Remove this line when you implement this function
            
        # Question B.2, here compute R1tF[i]
        pass  # TO DO... Remove this line when you implement this function

    # Question B.1, here compute the sample proportion p_Hmax, and the 95% confidence interval [lo_p_Hmax, hi_p_Hmax]
    pass  # TO DO... Remove this line when you implement this function
        
    print('p_Hmax = ' + format(p_Hmax,'.3e'))
    print('95% confidence interval = [' + format(lo_p_Hmax,'.3e') + ', ' + format(hi_p_Hmax,'.3e') + ']')
       
    # Question B.2, here compute the sample mean avg_R1tF, and the 95% confidence interval [lo_avg_R1tF, hi_avg_R1tF]
    pass  # TO DO... Remove this line when you implement this function

    print('avg_R1tF = ' + format(avg_R1tF,'.3e'))
    print('95% confidence interval = [' + format(lo_avg_R1tF,'.3e') + ', ' + format(hi_avg_R1tF,'.3e') + ']')
                
    plt.close('all')
    fig = plt.figure()
    Hbins = np.linspace(0.0,6e3,30)
    plt.hist(Hmax,bins=Hbins)
    plt.xlabel('Hmax')
    plt.ylabel('count')
    
    fig = plt.figure()
    plt.hist(R1tF,20)
    plt.xlabel('R1tF')
    plt.ylabel('count')
    plt.show()
    
    
    

    
    
if __name__ == "__main__":

    
# Comment out which test not to run
    
    plot_single_run()
    # stats_many_runs(nruns=1000)

    plt.show()
       

    
