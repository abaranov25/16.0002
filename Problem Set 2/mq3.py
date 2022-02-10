#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:16:22 2021

@author: allenbaranov
"""

import numpy as np
rng = np.random.default_rng()
N = 50000
X1 = rng.uniform(0.2, 0.8, int(N/2))
X2 = rng.uniform(0.0, 1.0, int(N/2))
X3 = np.concatenate([X1, X2])
Xr = np.random.permutation(X3) # This command re-orders the entries of X3 in a random manner

count = 0

for i in range(N):
    if (Xr[i] > .2): 
        count = count + 1
        
np.sum(Xr) / N
print("Part 1: {:d}".format(count))
print("Part 2: {:.3e}".format(count / N))
print("Part 3: {:.3e}".format(np.sum(Xr) / N))