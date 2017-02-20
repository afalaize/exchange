#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 12:43:23 2016

@author: Falaize
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


LAYOUT = 'circular'
ITERATIONS = 500


def node_color(node):
    return 'gainsboro'


def draw_nodes(graph, ax=None, layout=None, colors=None, root=None):
    if ax is None:
        ax = plt.axes(frameon=False)
    if colors is None:
        colors = [node_color(node) for node in graph.nodes()]
    if root is not None:
        graph.positions = hierarchy_pos(graph, root)
    if not hasattr(graph, 'positions'):
        if layout is None:
            layout = LAYOUT
        else:
            assert layout in ('circular', 'spring')
        if layout == 'spring':
            graph.positions = nx.spring_layout(graph, iterations=ITERATIONS)
        elif layout == 'circular':
            graph.positions = nx.circular_layout(graph)
    nx.draw_networkx_nodes(graph, graph.positions, ax=ax,
                           node_size=800, node_color=colors, lw=2)
    nx.draw_networkx_labels(graph, graph.positions, ax=ax)


def midle(nodes):
    return (nodes[0]+nodes[1])/2.


def vector(nodes):
    return nodes[1] - nodes[0]


def direction(v):
    return v/length(v)


def length(v):
    return np.sqrt(np.dot(v, v))


def orthogonal(v):
    mat = np.array([[0., 1.],
                    [-1., 0.]])
    return np.dot(mat, v)


def angle(d):
    return np.degrees((np.arctan(d[1]/d[0]) + np.pi/2 % np.pi) - np.pi/2)


def type_colors(type_=None):
    if type_ == 'storage':
        colors = ('blue', 'lightskyblue')
    elif type_ == 'dissipative':
        colors = ('green', 'lightsage')
    elif type_ == 'port':
        colors = ('red', 'lightsalmon')
    else:
        colors = ('k', 'gainsboro')
    return colors


def draw_edge(edge, pos, ax, move=0., forward=True, colors_type=None,
              draw=True):
    nodes = [pos[n] for n in edge[:2]]
    if colors_type is None:
        colors = type_colors()
    else:
        colors = type_colors(colors_type)
    if draw:
        conectionstyle = 'arc3, rad={0}'.format(move)
        patch = mpatches.FancyArrowPatch(*nodes,
                                         connectionstyle=conectionstyle,
                                         arrowstyle='wedge',
                                         mutation_scale=20.0,
                                         lw=2,
                                         edgecolor=colors[0],
                                         facecolor=colors[1])
        ax.add_patch(patch)

    bbox_props = dict(boxstyle="round, pad=0.3", fc=colors[1], ec=colors[0],
                      lw=2, alpha=1)
    v = vector(nodes)
    m = midle(nodes)
    d = direction(v)
    d_o = orthogonal(d)

    P0, P2 = nodes
    P1 = m + move*length(v)*d_o

    m1 = midle([P0, P1])
    m2 = midle([P1, P2])
    text_pos = midle([m1, m2])
    ax.text(list(text_pos)[0],
            list(text_pos)[1],
            edge[-1]['label'],
            ha="center", va="center",
            rotation=angle(d),
            size=12,
            bbox=bbox_props)


def moves(nodes, nedges, pos):
    MAX_ANGLE = np.pi/2
    m = list()
    if bool(nedges % 2):
        m.append(0.)
        nedges -= 1
    [m.append(((i+1)/(nedges+1))*MAX_ANGLE) for i in range(nedges//2)]
    [m.append(-((i+1)/(nedges+1))*MAX_ANGLE) for i in range(nedges//2)]
    return m


def draw_edges(graph, ax):
    for edge in graph.edges(data=True):
        draw_edge(edge, graph.positions, ax)


def plot(graph, filename=None, ax=None, layout=None, root=None):
    """
    plot of a PHSGraph
    """
    if ax is None:
        fig = plt.figure()
        ax = plt.axes(frameon=False)
    draw_nodes(graph, ax=ax, root=root)
    draw_edges(graph, ax)
    plt.tight_layout()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    plt.show()
    if filename is not None:
        if not filename[-4:] == '.pdf':
            filename += '.pdf'
        fig.savefig(filename)


def hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5 ):
    '''If there is a cycle that is reachable from root, then result will not be a hierarchy.

       G: the graph
       root: the root node of current branch
       width: horizontal space allocated for this branch - avoids overlap with other branches
       vert_gap: gap between levels of hierarchy
       vert_loc: vertical location of root
       xcenter: horizontal location of root
    '''

    def h_recur(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5,
                  pos = None, parent = None, parsed = [] ):
        if(root not in parsed):
            parsed.append(root)
            if pos == None:
                pos = {root:(xcenter,vert_loc)}
            else:
                pos[root] = (xcenter, vert_loc)
            neighbors = G.neighbors(root)
            if parent != None:
                try:
                    neighbors.remove(parent)
                except:
                    pass
            if len(neighbors)!=0:
                dx = width/len(neighbors)
                nextx = xcenter - width/2 - dx/2
                for neighbor in neighbors:
                    nextx += dx
                    pos = h_recur(G,neighbor, width = dx, vert_gap = vert_gap,
                                        vert_loc = vert_loc-vert_gap, xcenter=nextx, pos=pos,
                                        parent = root, parsed = parsed)
        return pos

    return h_recur(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5)