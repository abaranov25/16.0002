#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 13:26:02 2022

@author: allenbaranov
"""

t = [1,2,3]

s = t[::-1]

for element in zip(t, s):
    print(element)
    
s.append(t)
t.append(s)
    
for element in zip(t, s):
    print(element)
    

print(s)