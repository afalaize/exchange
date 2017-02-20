# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:07:22 2017

@author: afalaize
"""



from __future__ import absolute_import


import os
import numpy as np
from xml.etree import cElementTree as ElementTree
from decimal import Decimal


def array2text(data, dtype=None):
    if dtype is None:
        dtype = Decimal
    string = ''
    Nx, Nc = data.shape
    for i in range(Nx):
        formated_vector = map(dtype, data[i, :].tolist())
        string += ' '.join(['{:f}',]*3).format(*formated_vector) + '\n'
    return string
    

def getTemplate(filename):
    # Get file content
    vtufile = open(filename, mode='r')
    tree = ElementTree.parse(vtufile)
    vtufile.close()
    VTKFile = tree.getroot()
    UnstructuredGrid = VTKFile[0]
    Piece = UnstructuredGrid[0]
    for tag in ['CellData', 'UserData']:
        tags = [el.tag for el in Piece]
        index = tags.index(tag)
        del Piece[index]
    tags = [el.tag for el in Piece]
    PointData = Piece[tags.index('PointData')]
    names = [e.attrib['Name'] for e in PointData]
    for name in names:
        current_names = [e.attrib['Name'] for e in PointData]
        index = current_names.index(name)
        if name == 'Vitesse(m/s)':
            print(PointData[index].attrib)
        del PointData[index]
    return tree
    

ATTRIBUTES = {'NumberOfComponents': '3', 
              'type': 'Float32', 
              'Name': 'Vitesse(m/s)',
              'format': 'ascii'}    
    

def setData(tree, data, label='PodBasisElement', dtype=None):
    VTKFile = tree.getroot()
    UnstructuredGrid = VTKFile[0]
    Piece = UnstructuredGrid[0]
    tags = [el.tag for el in Piece]
    PointData = Piece[tags.index('PointData')]
    Nx, Nb, Nc = data.shape
    for b in range(Nb):
        print('Writing pod basis element {}/{}'.format(b+1, Nb))
        d = data[:, b, :]
        subel = ElementTree.SubElement(PointData, 'DataArray')
        attrib = ATTRIBUTES.copy()
        attrib.update({'Name': label+str(b+1)})
        for key in attrib.keys():
            subel.set(key, attrib[key])
        subel.text = array2text(d, dtype=dtype)
    return tree
    
    
if __name__ == '__main__':
    filename = 'F:\\TESTS_THOST\\170113_4_alimentations_sans_turbulence\\Re10000\\Results\\ThostA_000030.vtu'
    piece = getTemplate(filename)
    