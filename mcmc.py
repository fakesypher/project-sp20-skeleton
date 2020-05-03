import networkx as nx
import numpy as np
from mwrc_approx import mwrc_approx
from utils import is_valid_network, average_pairwise_distance
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



def metropolis_hastings(G, proposal_func, init_func, score_func,
                        num_iters, step=10):
    """
    Runs the metropolis-hastings algorithm for
    num_iters iterations, using proposal_func
    to generate samples and scorer to assign
    probability scores to samples.
      
    proposal_func -- function that proposes
        candidate state; takes in current state as
        argument and returns candidate state
    init_func -- function that proposes starting
        state; takes no arguments and returns a
        sample state
    score_func -- function that calculates f(y)/f(x)
        * g(y,x)/g(x,y); takes in two state samples
        (the current sample x then the candidate y).
    
    Returns a sequence of every step-th sample. You 
    should only sample on upon acceptance of a new
    proposal. Do not keep sampling the current state.
    
    Note the total number of samples will NOT be
    equal to num_iters. num_iters is the total number
    of proposals we generate.
    """
    
    current_state = init_func(G)
    sequence = [current_state]
    
    for i in range(num_iters):
        proposal = proposal_func(current_state, G)
        score = score_func(current_state, proposal)
        accept = np.random.uniform(0, 1)
        if min(1, score) >= accept:
            current_state = proposal
            if i % step == 0: 
                sequence.append(current_state)
            
    return sequence
    
def starting_state(G):
    return G

def sample_candidate(sample, G):
    num_nodes = G.number_of_nodes()
    node_in = list(sample.nodes())
    done = False
    tried = num_nodes
    node_tried = list(G.nodes())
    while not done:
        F = sample.copy()
        switch_try = np.random.randint(0, tried)
        tried = tried - 1
        switch = node_tried[switch_try]
        node_tried.remove(switch)

        if switch in node_in:
            F.remove_node(switch)
        else:
            F.add_node(switch)
            for neighbor in nx.all_neighbors(G, switch):
                if neighbor in node_in:
                    F.add_edge(switch, neighbor, weight=G[switch][neighbor]['weight'])
        if nx.is_dominating_set(G, F.nodes()) and nx.is_connected(F):
            done = True
    
    return F

def make_acceptance_scorer():
    def scorer(current_sample, candidate):
        def prob_f(graph):
            T = mwrc_approx(graph)
            return np.log(average_pairwise_distance(T))
        fy = prob_f(candidate)
        fx = prob_f(current_sample)
        return np.exp(fy - fx)
    return scorer
scorer = make_acceptance_scorer()

def mcmc(G, num_iter=10000):
    samples = metropolis_hastings(G, sample_candidate, starting_state, scorer, num_iter)
    return mwrc_approx(samples[-1]), samples
