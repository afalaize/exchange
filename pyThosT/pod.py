# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:15:13 2017

@author: afalaize
"""

from __future__ import division, absolute_import, print_function
import numpy as np
from tools import tensor2vector, sortIndices


class POD(object):

    def __init__(self, data):
        # store the 3-dimensionnal array for the data
        self.data_3d = data
        # number of space nodes (points), time steps and spatial components
        self.nx, self.nt, self.nc = self.data_3d.shape
        # store the 2-dimensionnal array for the data
        self.data_2d =tensor2vector(data)

    def splitMeanFluct(self, set_mean_to_zero=False):
        """
        Compute the mean and fluctuating fields (self.data_mean and
        self.data_fluc, respectively).
        
        The mean can be forced to 0 by passing the option set_mean_to_zero=True
        """
        mean = np.mean(self.data_2d, axis=1)
        if set_mean_to_zero:
            mean = 0*mean
        self.data_mean = mean
        self.data_fluc = np.zeros((self.nx*self.nc, self.nt))
        for i, col in enumerate(self.data_2d.T):
            self.data_fluc[:, i] = col - mean

    def computePODBasis(self, weighting_matrix="Id"):
        """
        nuild the POD basis from the snapshots (self.basis)
        """
        if weighting_matrix=="Id":
            W = 1.  # np.eye(self.Nx)
        else:
            assert False, 'Mass matrix is not defined'

        # Temporal correlation matrix
        self.C = np.dot(self.data_fluc.T, np.dot(W, self.data_fluc))

        # Temporal correlation matrix eigen decomposition
        eigen_vals, eigen_vecs = np.linalg.eig(self.C)

        # Remove the imaginary part (which should be numerically close to zero)
        eigen_vals = np.real(eigen_vals)
        eigen_vecs = np.real(eigen_vecs)
        
        # sort by decreasing eigen values
        indices = sortIndices(eigen_vals)
        
        self.C_eigen_vals = [eigen_vals[n] for n in indices]
        self.C_eigen_vecs = eigen_vecs[:, np.array(indices)]

        # Define POD basis
        self.basis = np.dot(self.data_fluc, self.C_eigen_vecs)
        self.npod = self.basis.shape[1]
    
    def computeModesEnergy(self):
        self.modes_energy = list()
        for i, val in enumerate(self.C_eigen_vals):
            mode_energy = sum(self.C_eigen_vals[:i+1])/sum(self.C_eigen_vals)
            self.modes_energy.append(mode_energy)

    def truncatePODBasis(self, threshold=1e-6):
        self.computeModesEnergy()
        self.threshold = threshold
        self.npod = [me >= 1-self.threshold for me in self.modes_energy].index(True)+1
        self.basis = self.basis[:, :self.npod]

    def checkPODBasisIsOrthonormal(self):
        M = np.abs(np.dot(self.basis.T, self.basis))
        for i in range(M.shape[0]):
            M[i, i] = 0
        print("val max out of diag from np.dot(basis.T, basis):{}".format(M.max()))

    def normalizePODBasis(self):
        for i in range(self.npod):
            scalarprod = np.dot(self.basis[:,i], self.basis[:,i])
            self.basis[:,i] = self.basis[:,i]/np.sqrt(scalarprod)
            
    def temporalCoefficients(self, Nmodes=None):
        if Nmodes is None:
            Nmodes = self.npod
        return np.dot(self.data_fluc.T, self.basis[:, :Nmodes])
    
    def reconstructData(self, Nmodes=None):
        if Nmodes is None:
            Nmodes = self.npod
        coeffs = self.temporalCoefficients(Nmodes=Nmodes)
        self.reconstruction = np.zeros(self.data_2d.shape)
        for t in range(coeffs.shape[0]):
            self.reconstruction[:, t] = self.data_mean + \
                np.dot(self.basis[:, :Nmodes], coeffs[t, :].T)
        
    def reconstructionError(self, Nmodes=None):
        if Nmodes is None:
            Nmodes = self.Npod
        self.reconstructData(Nmodes=Nmodes)
        
        def norm(data_2d):
            return np.mean((np.array(data_2d)*np.array(data_2d)).sum(axis=0))

        return norm(self.data_2d-self.reconstruction)
