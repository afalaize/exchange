# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 16:47:56 2017

@author: afalaize

- lecture de fichier xml au format ".vtu"
- les fichiers sont générés par ThosT
- le dossier "PATH" est un dossier "Results" de Thost et doit contenir 
    (i) les fichiers .vtu du cas courant
    (ii) un fichier ThostA.pvd qui répertorie tous les .vtu du cas courant 
"""


from __future__ import absolute_import


import os
import numpy as np
from xml.etree import cElementTree as ElementTree
from data_structures import Data, TimeSerie
from parallelization import concurentmap

PATH = r'F:\TESTS_THOST\170113_4_alimentations_sans_turbulence\Re10000\Results'

PVDNAME = r'ThostA.pvd'

DATANAMES = {'vitesse': 'Vitesse(m/s)',
             'pression': 'Pression(Bar)'}


def getElementTree(filename):
    """
    get the root ElementTree from the .vtu file with name "filename". 
    """
    vtufile = open(filename, mode='r')
    tree = ElementTree.parse(vtufile)
    vtufile.close()
    return tree

    
def getDataFromElementTree(tree, dataname=None):
    """
    Extract the data in xml file "filename" with tag "dataname" and return 
    the associated numpy.array
    
    Usage
    ------
    data = getDataFromVtu(filename, dataname=None)
    
    Parameters
    ----------
    filename: raw str
        The absolute path to the .vtu file (e.g. r"C:\folder\subfolder\ThostA_000051.vtu")
        
    dataname: raw str
        The label of the data to extract from the file pointed by "filename"
        
    Return
    -------
    text: string
        The data with tag "dataname" in file "filename" as a multilines string.
    """
    
    # name of the data to extract from the elementTree
    if dataname is None:
        dataname = 'Vitesse(m/s)'
    
    # Init output
    text = None
    
    # list all elements in the xml tree with tag DataArray
    for elem in tree.getiterator('DataArray'):
        # if the name of data is dataname
        if elem.attrib.has_key('Name') and elem.attrib['Name'] == dataname:
            # update text output with the xml element text
            text = elem.text
        
    # raise an error if text is still None (no DataArray with Name='dataname')
    assert text is not None, 'Can not find tag \n{} in xml tree.'.format(dataname)
    
    return text


def getCoordinatesFromElementTree(tree):
    """
    Return the coordinates of mesh points from xml tree (ElementTree) as text.
    """
    return tree.getiterator('Points')[0][0].text

    
def text2array(text):
    """
    format the data from text (string) to a 2D numpy.array.
    """
    # init the data as a list
    data = list()
    
    # get the data from each text line and append the data list
    for line in text.splitlines()[1:-1]:
        data.append(np.fromstring(line, dtype=float, sep=' '))
        
    # transform the list of 1D array to a 2D array and return
    return np.array(data)


def pvd2ListOfFiles(path, pvdname=None):
    """
    Extract a list of .vtu file names from the file pvdname in path.
    """
    if pvdname is None:
        pvdname = PVDNAME
    filename = path + os.sep + pvdname
    pvdtree = ElementTree.parse(open(filename))
    
    listOfFiles = list()
    for elem in pvdtree.getiterator('DataSet'):
        listOfFiles.append( path + os.sep + elem.attrib['file'])
    return listOfFiles


def elementTreesGenerator(listOfFiles):
    """
    Build a generator that reads and yields an elementTree for eac .vtu file 
    associated with the .pvd file with name 'pvdname' in 'path'
    """
    for filename in listOfFiles:
        yield getElementTree(filename)
    

def readTimeSerie(path, pvdname=None, dataname=None, count=None):
    """
    read all .vtu files in 'path' associated with .pvd file with name 'pvdname'
    (also in 'path'), and return a list of Data objects
    """

    if pvdname is None:
        pvdname = PVDNAME

    listOfFiles = pvd2ListOfFiles(path, pvdname=None)
    
    if count is None:
        count = len(listOfFiles)

    # name of the data array to extract from each vtu file
    if dataname is None:
        dataname = 'Vitesse(m/s)'

    def getDataFromTree(mesh_text, data_text, i):
        """
        function that yields a data_structures.Data for the data with name 
        'dataname' in the (element) tree.
        """
        print('Loading data from file {}/{}'.format(i+1, count))
        return Data(text2array(mesh_text), text2array(data_text))
    
    def argsGenerator():
        # generator for the ElementTrees in all files of the time serie
        tree_generator = elementTreesGenerator(listOfFiles[:count])
        i = 0
        for tree in tree_generator:
            data_text = getDataFromElementTree(tree, dataname)
            mesh_text = getCoordinatesFromElementTree(tree)
            yield (mesh_text, data_text, i)
            i += 1
            
    generator = argsGenerator()
    timeserie = map(lambda arg: getDataFromTree(*arg), generator)  
    print('Load of timeserie {} done.'.format(pvdname))
    return TimeSerie(tuple(timeserie))
