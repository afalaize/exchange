# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 17:04:26 2017

@author: afalaize
"""
from __future__ import division, print_function, absolute_import

import numpy
from scipy.interpolate import LinearNDInterpolator


def nodesCoordinatesToRegularGrid(nodes, h):
    minmax = list()
    npoints = list()
    ncomp = len(nodes[0,:])
    for i in range(ncomp):
        minmax.append([min(nodes[:, i]), max(nodes[:, i])])
        n = round((minmax[-1][1]-minmax[-1][0])/h)
        npoints.append(n)
    regular_axes = list()
    all_h = list()
    for (mi, ma), n in zip(minmax, npoints):
        regular_axes.append(numpy.linspace(mi, ma, n))
        all_h.append(regular_axes[-1][1]-regular_axes[-1][0])
    grid = mix(regular_axes[0], mix(regular_axes[1], regular_axes[2]))
    return grid, all_h


def interpData(original_data, original_coords, interp_coords):
    assert len(original_data.shape) == 1
    interp_func = LinearNDInterpolator(original_coords, original_data)
    return interp_func(interp_coords)
    
    
def mix(a, b):
    shape_a = a.shape
    shape_b = b.shape
    if len(shape_b) == 1:
        b = b.reshape((shape_b[0], 1))
    B = numpy.vstack([b, ]*shape_a[0])
    A = numpy.array([[ea, ]*shape_b[0] for ea in a]).reshape((shape_a[0]*shape_b[0], 1))
    return numpy.hstack((A, B))
    
