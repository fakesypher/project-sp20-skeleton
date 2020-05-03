import networkx as nx
from parse import read_input_file, write_output_file
from utils import is_valid_network, average_pairwise_distance
import sys
import os
import matplotlib.pyplot as plt
from networkx.algorithms.approximation import min_weighted_dominating_set
import numpy as np
import itertools
from mwrc_approx import mwrc_approx
from mcmc import mcmc

def solve(G):
    """
    Args:
        G: networkx.Graph

    Returns:
        T: networkx.Graph
    """
    if G.number_of_nodes() <= 1:
        return G
    if G.number_of_nodes() == 2:
        G.remove_node(0)
        return G
    
    apd = []
    T = naive_solve(G)
    if is_valid_network(G, T):
        apd.append(T)
    else: print('naive_solve not connected')

    T = greedy_mds(G)
    if is_valid_network(G, T):
        apd.append(T)
    else: print('greedy_mds not valid network')

    T = idea1(G, 'prune')
    if is_valid_network(G, T):
        apd.append(T)
    else: print('idea1 not valid network') 

    T = idea2(G)
    if is_valid_network(G, T):
        apd.append(T)
    else: print('idea2 not valid network')

    T = idea3(G)
    if is_valid_network(G, T):
        apd.append(T)
    else: print('idea3 not valid network')
    


    return min(apd, key=average_pairwise_distance)


    
def idea1(G, dom_set):
    """
    Construct new complete graph containing the sp between every node
    """
    all_sp_length = dict(nx.all_pairs_shortest_path_length(G))
    all_sp_path = dict(nx.all_pairs_shortest_path(G))
    complete = nx.complete_graph(G.number_of_nodes())
    for i in all_sp_length:
        for j in all_sp_length[0]:
            if i != j:
                sp_length = all_sp_length[i][j]
                complete[i][j]['weight'] = sp_length

    min_tree = mwrc_approx(complete)
    
    used = set()
    for edge in min_tree.edges():
        path = all_sp_path[edge[0]][edge[1]]
        prev = None
        curr = path[0]
        for node in path[1:]:
            prev = curr
            curr = node
            used.add((prev, curr))

    pruned_tree = G.edge_subgraph(list(used))

    if dom_set == 'prune':
        return _prune_degree_1(pruned_tree)

def idea2(G):
    """
    Use betweenness centrality
        sum of a-p-s-t that pass through v
    
    Concept: The higher the b-c is for a vertex, the more weight it contributes
    Idea: Greedily remove vertices that have high b-c but retain the connectedness of the graph
    Drawback: May leave unconnected - then still need to remove vertices
    """
    n = G.number_of_nodes()
    # True if the node is in the Tree or touching a node in the tree
    satisfied = [False for i in range(n)] 
    used = [False for i in range(n)]

    bc = nx.betweenness_centrality(G, n)
    # Nodes sorted by their betweenness centrality
    list_bc = [[k,v] for k, v in bc.items()]
    bc_order = list(sorted(list_bc, key=lambda x: x[1], reverse=True))
    
    H = G.copy()
    curr = 0
    no_used = 1
    while not all(satisfied):
        curr_node = bc_order[curr]
        neighbors = list(nx.all_neighbors(G, curr_node[0]))
        neighbor_used = 0
        for neighbor in neighbors:
            if used[neighbor]:
                neighbor_used += 1

        
        used[curr_node[0]] = True
        bc_order.remove(curr_node)
        satisfied[curr_node[0]] = True
        for neighbor in neighbors:
            satisfied[neighbor] = True
        curr -= 1
        

        curr += 1
        if curr >= len(bc_order):
            curr = 0
            no_used += 1
    
    H = G.copy()
    for index, val in enumerate(used):
        if not val:
            F = H.copy()
            F.remove_node(index)
            if nx.is_connected(F):
                H.remove_node(index)

    assert nx.is_connected(H), "H is not connected"

    return mwrc_approx(H)

def idea3(G):
    T, samples = mcmc(G, 100000)
    return T

    
def _prune_degree_1(T):
    assert nx.is_connected(T), 'tree must be connected'
    cop = T.copy()
    mod = T.copy()
    for node, degree in cop.degree():
        if degree == 1:
            mod.remove_node(node)
    return mod


def greedy_mds(G):
    mds = min_weighted_dominating_set(G)

    final = G.copy()
    for node in G.copy():
        if node not in mds:
            final.remove_node(node)
    
    if not nx.is_connected(final):
        return mwrc_approx(G)
    return mwrc_approx(final)


def naive_solve(G):
    """
    Find minimum spanning tree and prune off degree 1 nodes
    """
    T = mwrc_approx(G)
    if nx.is_connected(T):
        return _prune_degree_1(T)
    else:
        return _prune_degree_1(nx.minimum_spanning_tree(G))
    


# Here's an example of how to run your solver.

# Usage: python3 solver.py test.in

if __name__ == '__main__':
    path = './inputs/'
    for filename in os.listdir(path):
        if filename.endswith(".in"):
            G = read_input_file(path + filename)
            T = solve(G)
            assert is_valid_network(G, T)
            print(filename)
            print("Average  pairwise distance: {}".format(average_pairwise_distance(T)))
            write_output_file(T, 'outputs/{}.out'.format(filename.split('.')[0]))
# if __name__ == '__main__':
#     filename = 'small-1.in'
#     path = './inputs/'
#     G = read_input_file(path + filename)
#     T = solve(G)
#     plt.figure(1)
#     nx.draw(G, with_labels=True)
#     plt.figure(2)
#     nx.draw(T, with_labels=True)
#     plt.show()
#     print(min_weighted_dominating_set(G))
#     assert is_valid_network(G, T)
#     print("Average  pairwise distance: {}".format(average_pairwise_distance(T)))
#     write_output_file(T, 'test/{}.out'.format(filename.split('.')[0]))
# if __name__ == '__main__':
#     path = './inputs/'
#     for filename in os.listdir(path):
#         if filename.endswith(".in") and filename.startswith("small"):
#             G = read_input_file(path + filename)
#             T = solve(G)
#             assert is_valid_network(G, T)
#             print("Average  pairwise distance: {}".format(average_pairwise_distance(T)))
#             write_output_file(T, 'outputs/{}.out'.format(filename.split('.')[0]))
