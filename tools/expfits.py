# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:44:02 2022

@author: johanssm
"""

# Functions related to exponential distribution calculations.
 
import numpy as np
    

def expfit2(x, F):

# Exponential probability distribution fit, P=[1,sigma,x0,err]
# err = mean square fitting error
# (The parameter P[0] is always 1, used for compatibility with 
# the Weibull fit functions)
    xsolve = np.stack((x, -np.ones(len(x))), axis = 1)
    fsolve = -np.log(1 - F)
    Ps = np.linalg.lstsq(xsolve, fsolve, rcond = None)[0]
    s = 1 / Ps[0]
    x0 = Ps[1] / Ps[0]
    P = [1, s, x0]
    
    P.append(np.sum((cum_expdistr(P, x) - F) **2) / len(x))

    return P


def cum_expdistr(P, x):

# Cumulative exponential distribution

#    n = P[0]
    s = P[1]
    x0 = P[2]

    F = (1 - np.exp(-(x - x0) / s))

    return F
    
