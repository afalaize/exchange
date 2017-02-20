#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 14:05:15 2017

@author: Falaize
"""

from __future__ import division, absolute_import, print_function

import numpy as np
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl

from grid_tools import buildGrid

# Tweak how images are plotted with imshow
mpl.rcParams['image.interpolation'] = 'none' # no interpolation
mpl.rcParams['image.origin'] = 'lower' # origin at lower left corner
mpl.rcParams['image.cmap'] = 'RdBu_r'


class Data:
    """
    This is a structure for data with ncd components defined on a mesh in ncx
    dimensions with nx points.

    Parameters
    -----------

    mesh: 2d numpy.array with shape (nx, ncx)
    data: 2d numpy.array with shape (nx, ncd)

    Attributes
    -----------
    
    
    """
    def __init__(self, mesh, data):
        self.mesh = mesh
        self.nxm, self.ncx = mesh.shape

        nx, self.ncd = data.shape
        assert self.nxm == nx, 'number of mesh points {} is not equal to number\
 of data points {}.'.format(self.nxm, nx)
        self.data = data
        self.interpolator = LinearNDInterpolator(mesh, data)

    def interpToGrid(self, h=None, grid=None):
        self.grid_h = list()
        if h is None:
            assert grid is not None, "Either 'h' or 'grid' parameter should be\
 specified."
            self.grid_nd = grid

        else:   
            minmax = self.getMeshMinMax()
            self.grid_nd = self.buildGrid(minmax, h=h)

        self.grid_shape = self.grid_nd.shape
        self.grid = self.grid_nd.reshape((self.grid_shape[0],
                                          np.prod(self.grid_shape[1:]))).T
        for i, xi in enumerate(self.grid.T):
            self.grid_h.append(max(np.diff(xi)))

        self.data_grid = self.interpolator(self.grid)
        self.data_grid_nd = self.data_grid.T.reshape(([self.ncd, ] +
                                                      list(self.grid_shape[1:])
                                                      ))

    def computeGradient(self):
        gradient = [np.gradient(d, *self.grid_h) for d in self.data_grid_nd]
        self.gradient_nd = np.array(gradient)
        self.gradient = self.gradient_nd.reshape((self.ncd,
                                                  self.ncx,
                                                  np.prod(self.grid_shape[1:])
                                                  )).T
        self.gradient = np.swapaxes(self.gradient, 1, 2)

    def getMeshMinMax(self):
        minmax = list()
        for i, xi_mesh in enumerate(self.mesh.T):
            ximin, ximax = min(xi_mesh), max(xi_mesh)
            minmax.append((ximin, ximax))        
        return minmax
        
    @staticmethod
    def buildGrid(minmax, h=1.):
        return buildGrid(minmax, h=h)


# %%        
class TimeSerie:
    def __init__(self, dataList):
        self.serie = dataList
        self.nt = len(dataList)
        
    def interpOnGrid(self, grid):
        for i, d in enumerate(self.serie):
            print('interpolate data in TimeSerie: {}/{}'.format(i+1, self.nt))
            d.interpToGrid(grid=grid)

    def computeGradients(self):
        for i, d in enumerate(self.serie):
            print('compute gradient of data in TimeSerie: {}/{}'.format(i+1, self.nt))
            d.computeGradient()
            
    def formData3d(self, mode='grid'):
        assert mode in ('grid', 'mesh'), 'Mode {} not understood.'.format(mode)
        
        if mode == 'grid':
            attr = 'data_grid'
        elif mode == 'mesh':
            attr = 'data'
        data = list()
        for i, d in enumerate(self.serie):
            data.append(getattr(d, attr))
        self.data3D = np.swapaxes(np.array(data), 0, 1)
        
# %%
def test2D():
    nx = 10000

    nx1 = int(nx**(1/2.))
    x1 = np.linspace(0, 1, nx1)

    nx2 = int((nx/nx1)**(1.))
    x2 = np.linspace(0, 1, nx2)

    X = np.array(np.meshgrid(x1, x2))

    nx = nx1*nx2
    mesh = X.reshape((2, nx)).T


    def f_2d(x, y):
        '''a function with 2D input to interpolate on [0,1]'''
        return np.exp(-x)*np.cos(x*2*pi)*np.sin(y*2*pi)

    f_2d_grid = f_2d(x1.reshape(-1, 1), x2)
    f_2d_grid = f_2d_grid.reshape((1, nx)).T

    print(f_2d_grid.shape)

    data = Data(mesh, f_2d_grid)
    data.interpToGrid(h=0.001)

    data.computeGradient()
    g = data.gradient

    plt.figure()
    plt.imshow(data.data.reshape((nx1, nx2)).T)
    plt.figure()
    plt.imshow(data.data_grid_nd[0].T)

    plt.figure()
    plt.imshow(g[:, 0, 1].reshape(data.grid_shape[1:]).T)
    plt.figure()
    plt.imshow(data.gradient_nd[0][1].T)
    
    return data


# %%
def test3D():
    nx = 10000

    nx1 = int(nx**(1/3.))
    x1 = np.linspace(0, 1, nx1)

    nx2 = int((nx/nx1)**(1/2.))
    x2 = np.linspace(0, 1, nx2)

    nx3 = int((nx/nx1/nx2))
    x3 = np.linspace(0, 1, nx3)

    X = np.array(np.meshgrid(x1, x2, x3))

    nx = nx1*nx2*nx3
    mesh = X.reshape((3, nx)).T

    def f_3d(x, y, z):
        '''a function with 3D input to interpolate on [0,1]'''
        return np.sin(x*2*pi)*np.sin(y*2*pi)*np.sin(z*2*pi)

    f_3d_grid = f_3d(x1.reshape(-1, 1, 1), x2.reshape(1, -1, 1), x3)
    f_3d = f_3d_grid.reshape((nx, 1))

    data = Data(mesh, f_3d)

    data.interpToGrid(h=0.01)
    return data


# %%

