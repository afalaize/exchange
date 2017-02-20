# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 16:05:42 2017

@author: afalaize
"""

def indent(string):
    "Prepend each line of the string with four empty spaces."
    return "\n".join(['    ' + el for el in string.split('\n')])