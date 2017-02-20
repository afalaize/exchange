# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 13:51:20 2017

@author: afalaize
"""


ELEMENTS = {'P0',
            'P0C',
            'P1'}

QUANTITITES = {'Scalaire',
               'Vecteur',
               'Tenseur',
               'Tenseur_Sym'}

class Champ:
    def __init__(self, **kwargs):
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        