################################################################################
# 16.0002 / 18.0002 Fall 2021
# Problem Set 3: myrootfinding
# Name: Allen Baranov
# Collaborators:

import numpy as np

def bisection(f,a,b,n):

    """
    Runs the bisection method for finding an approximation to a root of f
    
    Args:
	f (function object): scalar function
        a, b (float): numbers such that f(a)*f(b) < 0
		(i.e., f(a) and f(b) have opposite signs)
        n (int): number of iterations of the bisection method

    Returns:
        x: approximation to a root of f, i.e., such that f(x) = 0

    """
        
    # Here program the bisection method
    
    x = 0
    a1 = a
    b1 = b
    

    for i in range(n):
        x = (a1+b1)/2
        if (f(x)*f(a1) > 0):
            a1 = x
        else:
            b1 = x
    
    return x



