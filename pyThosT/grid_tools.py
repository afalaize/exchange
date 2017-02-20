# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 10:48:10 2017

@author: afalaize
"""

from __future__ import division, print_function, absolute_import
import numpy as np


def buildGrid(minmax, h=1.):
    """
    Return an N-dimensional regular grid.
    
    Parameters
    -----------
    minmax: iterable 
        N tuples of 2 floats: minmax = [(x1_min, x1_max), ..., (xN_min, xN_max)] 
        with (xi_min, xi_max) the grid limits over the i-th axis.
    
    h: float or iterable
        Grid spacing(s). If h is a float, the grid spacing is the same over every 
        axis. If h is an iterable with length N, h[i] is the grid spacing over
        axis i.
        
    Return
    -------
    grid: array
        An N-dimensional regular grid.
             
    """
    x_grid = list()
    for i, xminmax in enumerate(minmax):
        ximin, ximax = xminmax
        hi = h[i] if isinstance(h, (tuple, list)) else h
        nxi = int((ximax-ximin)/hi)
        xi_grid = np.linspace(ximin, ximax, nxi)
        x_grid.append(xi_grid)
    return np.array(np.meshgrid(*x_grid))
