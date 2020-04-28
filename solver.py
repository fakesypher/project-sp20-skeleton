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
    apd = []
    T = naive_solve(G)
    apd.append(T)
    #T = minimum_routing_cost_approx(G)
    #apd.append(T)
    T = greedy_mds(G)
    apd.append(T)
    
    return min(apd, key=average_pairwise_distance)

def minimum_routing_cost_approx(G):
    """
    Returns the approximate minimum routing cost tree
    """
    return ...

def greedy_mds(G):
    nodes = [False for node in G.nodes()]
    deg_sort = sorted(G.degree(), reverse=True, key=lambda x: x[0])
    T = nx.Graph()
    while not all(nodes):
        top = deg_sort.pop()
        T.add_node(top[0])
        for node in nx.all_neighbors(G, top[0]):
            nodes[node] = True
    T_nodes = T.nodes()
    final = G.copy()
    for node in G.copy():
        if node not in T_nodes:
            final.remove_node(node)

    if not nx.is_connected(final):
        return nx.minimum_spanning_tree(G)
    return nx.minimum_spanning_tree(final)


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
