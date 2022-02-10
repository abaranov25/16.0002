################################################################################
# 16.0002 / 18.0002 Fall 2021
# Problem Set 2: mysolve.py
# Name: Allen Baranov
# Collaborators:

"""
Solve Ku = f by Gaussian elimination (no pivoting)

"""
import numpy as np

def mysolve(K, f):
    """    

    Args:
        K (numpy ndarray): square matrix
        f (numpy ndarray): right-hand side vector
            
    Returns:
        numpy ndarray: u, solution to Ku=f

    """
    m,n = K.shape
    assert m == n, "Non-square matrix"
    
    u = np.zeros(n)
    
    # Extended matrix with the right-hand side as last column
    A = np.zeros((n,n+1))
    A[:,0:n] = K
    A[:,n] = f
    
    # Elimination of nonzero elements below the diagonal
    for i in range(n):
        assert A[i,i] != 0.0, "Zero pivot detected"
        
        for j in range(i+1,n):
            Lij = A[j,i] / A[i,i]
            A[j,:] = A[j,:] - Lij * A[i,:]
    
    # Back substitution
    u[n-1] = A[n-1,n]/A[n-1,n-1]
    
    for i in range(n-2,-1,-1):
        u[i] = A[i,n]
        
        for j in range(i+1,n):
            u[i] = u[i] -  A[i,j]*u[j]
            
        u[i] = u[i] / A[i,i]
        
    return u
                
        
            
def mysolve_vectorized(K, f):
    """    

    Args:
        K (numpy ndarray): square matrix
        f (numpy ndarray): right-hand side vector
            
    Returns:
        numpy ndarray: u, solution to Ku=f

    """
    m,n = K.shape
    assert m == n, "Non-square matrix"
    
    u = np.zeros(n)
    
    # Extended matrix with the right-hand side as last column
    A = np.zeros((n,n+1))
    A[:,0:n] = K
    A[:,n] = f
    
    # Elimination of nonzero elements below the diagonal
    for i in range(n):
        assert A[i,i] != 0.0, "Zero pivot detected"
        
        for j in range(i+1,n):
            Lij = A[j,i] / A[i,i]
        
            ##### REPLACE THE TWO LINES BELOW #####
            #            for k in range(n+1):
            #                 A[j,k] = A[j,k] - Lij * A[i,k]
            
            A[j] = A[j] - Lij * A[i]
                    
    # Back substitution
    u[n-1] = A[n-1,n]/A[n-1,n-1]
    
    for i in range(n-2,-1,-1):
    
        u[i] = A[i,n]

            ##### REPLACE THE TWO LINES BELOW #####
	    #        for j in range(i+1,n):
	    #            u[i] = u[i] -  A[i,j]*u[j]
        
        u[i] = u[i] - np.dot( A[i][i+1:n] , u[i+1:] )

                    
        u[i] = u[i] / A[i,i]
        
    return u    