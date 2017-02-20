# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:53:14 2017

@author: afalaize
"""

from __future__ import division, print_function
import os


# CONFIGURATION
EXTENSION = r'mtc'


def read(filename, extension=None):
    """
    Read data from 'filename' and return the corresponding string.
    If extension is not provided, use 'mtc' as file extension.
    If filename does not end with '.'+extension, this is appended.
    """
    if extension is None:
        extension = EXTENSION
    if not filename.endswith(r'.'+extension):
        filename += r'.'+extension
    string = str()
    file_ = open(filename, "r")
    with file_ as openfileobject:
        for line in openfileobject:
            string += line
    file_.close()
    return string


def get_files(path, extension=None):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    if extension is None:
        extension = EXTENSION

    file_paths = []  # List which will store all of the full filepaths.
    # Walk the tree.
    for root, directories, files in os.walk(path):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            if filepath.endswith(r'.'+extension):
                file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.
