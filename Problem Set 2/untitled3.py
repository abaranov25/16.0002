#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  7 21:13:18 2021

@author: allenbaranov
"""

import numpy as np
#import yfinance as yf
import matplotlib.pyplot as plt

A = [1,2,3,4,5,6,7,8,9,10,11,12]


################################################################################
# Confidence interval: 
# 
#   Input: A (numpy array), n
#
#   Output: u (float list) of lower limit and upper limit of 95$ confidence
################################################################################




def ci(A):
    m = np.mean(A)
    s = np.std(A)
    
    u = []
    u.append(m-2*s)
    u.append(m+2*s)
    return u




if (A[-1] > ci(A)[1]):
    print("Alert")
    
    
    
    
    
x = np.arange(200)
y = np.arange(100)

z = np.zeros((100,200))
for j in range(200):
    for i in range(100):
        z[i,j] = i+j

fig, axs = plt.subplots()
c = axs.contour(x,y,z)

axs.clabel(c)
axs.grid(True)
axs.axis('square')