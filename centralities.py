import numpy as np
import os 
from networkx import from_numpy_array, to_numpy_array
import datetime
import itertools
from operator import itemgetter 

#centrality measures
from networkx.algorithms.centrality import degree_centrality, eigenvector_centrality, closeness_centrality, betweenness_centrality, betweenness_centrality_subset, katz_centrality

def betweenness_c(G, residue_names_1, n=10):
    """
    Compute the betweenness centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the betweenness centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: P_coef[node]}, for each node is linked its betweenness centrality
    """
    bc = betweenness_centrality(G)
    #bc= betweenness_centrality_parallel(G)
    bc = {int (float (k)):v for k,v in bc.items()}
    dict_node_centrality = dict ()

    for i, cent in bc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_bc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by betweenness centrality".format(n))
    for d in sorted_bc[:n]:
        print(d)

    return dict_node_centrality

def eigenvector_c(G, residue_names_1, n=10):
    """
    Compute the eigenvector centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the eigenvector centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: eigenvector_centrality[node]}, for each node is linked its eigenvector centrality
    """
    ec = eigenvector_centrality(G, max_iter=10000)
    ec = {int (float (k)):v for k,v in ec.items()}
    dict_node_centrality = dict ()

    for i, cent in ec.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_ec = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by eigenvector centrality".format(n))
    for d in sorted_ec[:n]:
        print(d)

    return dict_node_centrality

def degree_c(G, residue_names_1, n=10):
    """
    Compute the degree centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the degree centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: degree_centrality[node]}, for each node is linked its degree centrality
    """
    dc = degree_centrality(G)
    dc = {int (float (k)):v for k,v in dc.items()}
    dict_node_centrality = dict ()

    for i, cent in dc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_dc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by degree centrality".format(n))
    for d in sorted_dc[:n]:
        print(d)

    return dict_node_centrality

def closeness_c(G, residue_names_1, n=10):
    """
    Compute the closeness centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the closeness centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: closeness_centrality[node]}, for each node is linked its closeness centrality
    """
    cc = closeness_centrality(G)
    cc = {int (float (k)):v for k,v in cc.items()}
    dict_node_centrality = dict ()

    for i, cent in cc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_cc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by closeness_centrality".format(n))
    for d in sorted_cc[:n]:
        print(d)

    return dict_node_centrality

def katz_c(G, residue_names_1, n=10):
    """
    Compute the katz centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the katz centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: closeness_centrality[node]}, for each node is linked its katz centrality
    """
    kc = katz_centrality(G, max_iter=10000)
    kc = {int (float (k)):v for k,v in kc.items()}
    dict_node_centrality = dict ()

    for i, cent in kc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_kc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by katz centrality".format(n))
    for d in sorted_kc[:n]:
        print(d)

    return dict_node_centrality

def participation_coefs(G, labels, residue_names_1):
    """
    Compute the participation coefficient of the nodes of the graph given a partition.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the closeness centrality.
        labels: np.array, extracted clusters/communities
        residue_names_1: np.array, list of the residues names of the protein.
    Returns:
        P: dict {node: P_coef[node]}, for each node is linked its participation coefficient
    """
    A = to_numpy_array(G)
    n = A.shape[0]
    P = dict()
    k_s = np.zeros((n))

    for i in range(n):
        k_i = np.sum(A[i,:])
        k_si = 0

        for j in range(n):
            if (i!=j):
                if ((labels[i] == labels[j]) and (A[i,j]!=0)):#se il nodo i e il nodo j sono dello stesso cluster e c'è un arco che li connette
                    k_si += A[i,j]

        k_s[i] = k_si
        P[residue_names_1[i]] = 1 - (k_s[i]/k_i)**2

    return P