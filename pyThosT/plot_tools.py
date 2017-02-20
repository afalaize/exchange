# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 16:48:27 2017

@author: afalaize
"""
import matplotlib.pyplot as plt
import os
import numpy as np


def plot_kinetic_energy(energy, PATH):
    plt.figure()
    plt.plot(energy)
    plt.title('energie cinetique $E(t_j)=\sum_i (u(x_i, t_j)^T\cdot u(x_i, t_j))/2$')
    plt.xlabel('$j$')
    plt.ylabel('$E(t_j)$')
    plt.savefig(PATH + os.sep + 'energie_cinetique.pdf')

    
def plot_kinetic_energy_variation(energy, PATH):
    plt.figure()
    plt.plot(np.diff(energy))
    plt.title("variation d'energie cinetique $\Delta E(t_j) = E(t_{j})-E(t_{j-1})$")
    plt.xlabel('$j$')
    plt.ylabel('$\Delta E(t_j)$')
    plt.savefig(PATH + os.sep + 'variation_energie_cinetique.pdf')    


def plot_all_reconstruction_errors(pod, PATH):
    reconstruction_errors = list()
    for n in range(0, pod.npod):
        print("Compute reconstruction error for {} modes (out of {})".format(n+1, pod.npod))
        reconstruction_errors.append(pod.reconstructionError(Nmodes=n+1))
    plt.figure()
    plt.plot(range(1 ,pod.npod), reconstruction_errors[:-1])
    plt.title("erreur de reconstruction (!= eq (2.64) de la these d'Erwan)")
    plt.xlabel('Nombre de modes retenus')
    plt.ylabel('Erreur')
    plt.savefig(PATH + os.sep + 'erreur_de_reconstruction.pdf')
    

def plot_pod_eigenvals(pod, PATH):
    plt.figure()
    pod.compute_POD_modes_enrgy()
    plt.plot(1-np.array(pod.POD_modes_energy))
    plt.title("Valeurs propres de l'operateur de correlation")
    plt.xlabel('indice')
    plt.ylabel('Valeur')
    plt.savefig(PATH + os.sep + 'valeurs_propre.pdf')
            
