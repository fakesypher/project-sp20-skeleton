import networkx as nx
from parse import read_input_file, write_output_file
from utils import is_valid_network, average_pairwise_distance
import sys
import os

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
    T = naive_solve(G)
    
    return T

def minimum_routing_cost_approx(G):
    """
    Returns the approximate minimum routing cost tree
    """



def naive_solve(G):
    """
    Find minimum spanning tree and prune off degree 1 nodes
    """
    T = nx.minimum_spanning_tree(G)
    cop = T.copy()
    for node, degree in cop.degree():
        if degree == 1:
            T.remove_node(node)
    return T


# Here's an example of how to run your solver.

# Usage: python3 solver.py test.in

if __name__ == '__main__':
    path = './inputs/'
    for filename in os.listdir(path):
        if filename.endswith(".in"):
            G = read_input_file(path + filename)
            T = solve(G)
            assert is_valid_network(G, T)
            print("Average  pairwise distance: {}".format(average_pairwise_distance(T)))
            write_output_file(T, 'outputs/{}.out'.format(filename.split('.')[0]))
