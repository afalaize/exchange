# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:31:15 2017

@author: afalaize
"""

from __future__ import division, print_function

from files_io import read, get_files
from blocs import (get_bloc_limits, get_blocs_recursively)
import networkx as nx
from plot_graphs import plot as nxplot
import copy


ATTRIBUTES = {'Type',
              'Data',
              'Dependance',
              'DependanceModifiable',
              'ModeleTerminaison'}


class MainClass:
    def __init__(self, path):
        self.all_files = get_files(path)
        self.all_mtc = list()
        self.graph = nx.DiGraph()
        for filename in self.all_files:
            mtc = Mtc()
            mtc.read_file(filename=filename)
            mtc.get_all_models()
            self.all_mtc.append(mtc)
            self.graph.add_edges_from(mtc.graph.edges(data=True))


class Mtc:
    def __init__(self):
        self.string = r''
        self.models_str = list()
        self.models = list()
        self.graph = nx.DiGraph()

    def read_file(self, filename):
        "Read the file pointed by filename and store the data in self.string"
        self.filename = filename
        self.string += read(filename)

    def get_all_models_string(self):
        continuation_flag = True
        string = self.string
        while continuation_flag:
            limits = get_bloc_limits(string)
            if not limits == -1:
                start, end = limits
                self.models_str.append(string[start:end])
                string = string[end:]
            else:
                continuation_flag = False

    def get_all_models(self):
        self.get_all_models_string()
        for string in self.models_str:
            model = Model(string)
            self.models.append(model)
            self.graph.add_edges_from(model.graph.edges(data=True))

    def names(self):
        names = list()
        for model in self.models:
            names.append(model.Name)
        return names

    def types(self):
        types = list()
        for model in self.models:
            types.append(model.Type)
        return types


class Model:
    def __init__(self, string):
        # the raw text taken from the .mtc file
        self.string = string
        # the raw text taken from the .mtc file
        self.dico = get_blocs_recursively(string)
        self.Name = self.dico.keys()[0]
        for name in ATTRIBUTES:
            setattr(self, name, None)
            if isinstance(self.dico[self.Name], list):
                for item in self.dico[self.Name]:
                    if item.keys()[0] == name:
                        attr = item[name]
                        setattr(self, name, attr)
        self.buildGraph()

    def buildGraph(self):
        self.graph = nx.DiGraph()
        edges = list()
        root = self.Name
        if root == 'Lanceur':
            edge = (root, self.dico[self.Name], {'label': 'START', 'color':'b'})
            edges.append(edge)
        if self.Dependance is not None:
            if isinstance(self.Dependance, dict):
                node = self.Dependance.values()[0]
                edge = (root, node, {'label': self.Dependance.keys()[0]})
                edges.append(edge)
            else:
                reverse_deps = copy.copy(self.Dependance)
                reverse_deps.reverse()
                for i, dependance in enumerate(reverse_deps):
                    N1 = root
                    N0 = dependance.values()[0]
                    if self.Type == 'ModeleDeModeles' and i > 0:
                        N1 = N0
                    c = 'r' if self.Type == 'ModeleDeModeles' else 'k'
                    edge = (N1, N0, {'label': dependance.keys()[0],
                                     'color': c})
                    edges.append(edge)
        self.graph.add_edges_from(edges)


if __name__ == '__main__':
    path = r'F:\TESTS_THOST\161215_Ahmed\CasTest1'
#    path = '/Users/Falaize/Documents/RECHERCHE/RSYNC_LASIE/DEV/MTC/161215_DocsAhmed/formation_cimlib/5_remaillage/3fois'
    klass = MainClass(path)
    g = klass.graph
    nxplot(g, root='Lanceur')
