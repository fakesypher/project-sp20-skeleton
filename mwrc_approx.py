import networkx as nx
from parse import read_input_file, write_output_file
from utils import is_valid_network, average_pairwise_distance
import sys
import os
import matplotlib.pyplot as plt
from networkx.algorithms.approximation import min_weighted_dominating_set
import numpy as np
import itertools
import math

def mwrc_approx(G):
    """
    Compos's algorithm
    https://doi.org/10.1016/j.comnet.2008.08.013
    """

    if G.number_of_edges() <= 1:
        return G

    C_1, C_2, C_3 = .2, .2, .6
    nodes = G.nodes()
    node_id = {node: i for i, node in enumerate(nodes)}
    id_node = {i: node for i, node in enumerate(nodes)}
    #
    # Compute Degree and Weights for each vertex
    #
    d = [0 for _ in nodes] # degree
    s = [0 for _ in nodes] # sum of adjacent edge
    m = [0 for _ in nodes] # max of adjacent edges
    total_weight = 0 # total_weight
    for edge in G.edges():
        i, j = node_id[edge[0]], node_id[edge[1]]
        
        w = G[edge[0]][edge[1]]['weight']
        d[i] += 1
        s[i] += w
        m[i] = max(m[i], w)
        d[j] += 1
        s[j] += w
        m[j] = max(m[j], w)
        total_weight += w

    w = [101 for _ in nodes] # weight
    cf = [100001 for _ in nodes] # estimated cost to f
    sp = [0 for _ in nodes] # spanning potential
    sp_max = 0 # max spanning potential
    f = 0
    p = [None for _ in nodes] # parents
    pd = [None for _ in nodes] # parent degree
    ps = [None for _ in nodes] # sum of parent neighbor weight


    #
    # Calculate Mean and StdDev and set C_4, C_5
    #
    mean = total_weight / G.number_of_edges()
    tot = 0
    for edge in G.edges():
        tot = tot + (G[edge[0]][edge[1]]['weight'] - mean)**2
    std = math.sqrt(tot / (G.number_of_edges() - 1))
    ratio = std / mean
    C_4, C_5 = 0, 0
    if ratio < .4 + .005 * (len(nodes) - 10):
        C_4, C_5 = 1, 1
    else:
        C_4, C_5 = .9, .1

    #
    # Select highest sp as initial vertex
    #
    for un_node in nodes:
        node = node_id[un_node]
        sp[node] = C_1 * d[node] + C_2 * d[node] / s[node] + C_3 / m[node]
        if sp[node] > sp_max:
            f = node
            sp_max = sp[node]
    
    w[f] = 0
    cf[f] = 0
    p[f], pd[f], ps[f] = f, 0, 1
    track = [None for _ in nodes]

    L = set([f]) # initialize set L
    T = set()
    track[f] = 1 # 1 means f is in L, 2 means f is in T
    num_spanned = 0
    wd = [1000001 for _ in nodes]
    jsp = [0 for _ in nodes]
    wd[f], jsp[f] = C_4 * w[f] + C_5 * cf[f], d[f] + pd[f] + (d[f] + pd[f])/(s[f] + ps[f])
    
    while num_spanned < len(nodes):
        
        wd_min, jsp_max = 100001, 0
        u = 0
        

        for node in L:
            if wd[node] < wd_min:
                S = set([node])
                wd_min = wd[node]
            elif wd[node] == wd_min:
                S.add(node)


        for node in S:
            if jsp[node] >= jsp_max:
                jsp_max = jsp[node]
                u = id_node[node]
        L.remove(node_id[u])

        for neighbor in nx.all_neighbors(G, u):
            v = node_id[neighbor]
            u_id = node_id[u]
            if track[v] == 2: # 2 means in T already
                continue
            uv_weight = G[u][neighbor]['weight']
            wd_new = C_4 * uv_weight + C_5 * (cf[u_id] + uv_weight)
            dv, du = d[v], d[u_id]
            jsp_new = dv + du + (dv + du) / (s[v] + s[u_id])
            if wd_new < wd[v]:
                wd[v] = wd_new
                jsp[v] = jsp_new
                p[v] = u_id
                cf[v] = cf[u_id] + G[u][neighbor]['weight']
                pd[v] = d[u_id]
                ps[v] = s[u_id]
            elif (wd_new == wd[v] and jsp_new >= jsp[v]):
                jsp[v] = jsp_new
                p[v] = u_id
                cf[v] = cf[u_id] + G[u][neighbor]['weight']
                pd[v] = d[u_id]
                ps[v] = s[u_id]

            if track[v] == None:
                L.add(v)
                track[v] = 1

        track[u_id] = 2
        num_spanned += 1
        par = id_node[p[u_id]]
        if par != u:
            T.add((u, par, G[u][par]['weight']))

    tree = nx.create_empty_copy(G)
    tree.add_weighted_edges_from(list(T))
    return tree



