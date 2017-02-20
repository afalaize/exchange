# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:17:03 2017

@author: afalaize
"""


import numpy as np


def produit_scalaire(v1, v2):
    """
    Inputs
    ------
    
    v1, v2: np.array with shape (nb nodes, nb components)
    
    Output
    -------
    
    result: float

    """
    
    def columns_generator(array):
        """
        Yield the columns of array
        """
        for i, c in enumerate(np.transpose(array)):
            yield c
            
    return sum(
               map(lambda u, v: np.dot(u, v), 
                   columns_generator(v1), columns_generator(v2)
                   )
               )
  
               
def trace(array):
    """
    Inputs
    ------
    
    array: np.array with shape (nb nodes, nb components, nb components)
        Order 2 tensor defined over a grid.
    
    Output
    -------
    
    result: 1d np.array with shape (nb nodes, ).

    """
    
    nx, nc, nc_bis = array.shape
    
    assert nc == nc_bis
    
    def tensor_generator():
        for _, e in enumerate(array):
            yield e
            
    return np.array(map(lambda t: t.trace(), tensor_generator()))
    