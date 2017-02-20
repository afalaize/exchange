# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 15:45:31 2017

@author: afalaize
"""

import os
import numpy as np

from pod import POD
from vtu_reader import readTimeSerie, pvd2ListOfFiles
from vtu_writer import getTemplate, setData
from plot_tools import plot_kinetic_energy, plot_kinetic_energy_variation
from tools import compute_kinetic_energy, vector2tensor
from grid_tools import buildGrid

CONFIG = {'path': r'F:\TESTS_THOST\170113_4_alimentations_avec_turbulence\Re50000\Results',
          'threshold': 1e-9,
          'count': None,
          'set_mean_to_zero': False,
          'outputfilename': 'POD.vtu',
          'h_interp': 1e-2
          }

          
# Recover the data from the timeserie (pvd)
TS = readTimeSerie(CONFIG['path'], count=CONFIG['count'])

#%%
minmax = TS.serie[0].getMeshMinMax()
grid = buildGrid(minmax, CONFIG['h_interp'])          

TS.interpOnGrid(grid)
TS.formData3d(mode='grid')

#%%

pod = POD(TS.data3D)
energy = compute_kinetic_energy(pod.data_3d)

pod.splitMeanFluct(set_mean_to_zero=CONFIG['set_mean_to_zero'])
pod.computePODBasis()

pod.truncatePODBasis(threshold=CONFIG['threshold'])
pod.normalizePODBasis()
pod.checkPODBasisIsOrthonormal()

#basis_3d = vector2tensor(pod.basis, 3)
#
#listOfFiles = pvd2ListOfFiles(CONFIG['path'])
#templateVtu = getTemplate(listOfFiles[-1])
#templateVtu = setData(templateVtu, basis_3d, dtype=float)
#
#outputfilename = CONFIG['path'] + os.sep + CONFIG['outputfilename']
#outputfile = open(outputfilename, 'w')
#templateVtu.write(outputfile)
#outputfile.close()
#
