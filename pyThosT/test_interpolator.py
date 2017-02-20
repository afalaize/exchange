# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 11:35:23 2017

@author: afalaize
"""

from interpolation import interpolator, NodesCoordinatesToRegularGrid
from vtu_reader import getElementTree, readNodesCoordinates

FILENAME = r'F:\TESTS_THOST\170113_4_alimentations_avec_turbulence\Re50000\Results\ThostA_000081.vtu'

coordinates = readNodesCoordinates(getElementTree(FILENAME))

grid, all_h = NodesCoordinatesToRegularGrid(coordinates, h=1e-2)


#%%
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.scatter(grid[0:-1,0], grid[0:-1,1], grid[0:-1,2], c='b', marker='.')