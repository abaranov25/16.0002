Problem Set 1: Applied IVP Implementation

In this problem set, I implemented an IVP (initial value problem) solver by creating a class and helper methods (both the Runge-Kutta 4 algorithm and 
Euler Method). An IVP problem is a differential equation that cannot be manually solved, so the idea is to take timesteps given an initial set of values 
and an equation defining the rate of change (ODE or PDE) as a function of the current state. 

This was applied to solve for the trajectory of a Martian lander and a hail particle, and plots were made to illustrate the calculated trajectories using
the matplotlib library. 
