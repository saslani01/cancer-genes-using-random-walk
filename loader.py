import numpy as np
from scipy.sparse.csgraph import connected_components

# lcc (largest connected component)
def protein_interaction_adj_matrix(path, lcc=False):
    protein_index = {}
    interactions = []

    with open(path, 'r') as f:
        # index the proteins and store interactions
        for line in f:
            p1, p2 = line.strip().split()
            if p1 not in protein_index:
                protein_index[p1] = len(protein_index)
            if p2 not in protein_index:
                protein_index[p2] = len(protein_index)
            interactions.append((p1, p2))

    # n * n zero matrix
    n = len(protein_index)
    adj_matrix = np.zeros((n, n))
    
    # fill the adj matrix 
    for p1, p2 in interactions:
        adj_matrix[protein_index[p1]][protein_index[p2]] = 1
        adj_matrix[protein_index[p2]][protein_index[p1]] = 1

    if lcc:
        labels = connected_components(adj_matrix, directed=False)[1]
        lcc_indices = np.where(labels == np.argmax(np.bincount(labels)))[0] # need [0] because the result is inside a tuple

        adj_matrix = adj_matrix[np.ix_(lcc_indices, lcc_indices)]
        lcc_protein_index = {name: i for name, i in protein_index.items() if i in lcc_indices} # filter
        protein_index = {name: i for i, name in enumerate(lcc_protein_index)} # reindex

        # checking if our lcc contains all the starting onco genes
        missing = [g for g in cancer_genes("./data/onco_genes.txt") if g not in protein_index]
        if missing:
            print(f"Onco Genes Missing in LCC: {missing}")

    return adj_matrix, protein_index

def cancer_genes(path):
    genes = []
    with open(path, 'r') as f:
        for line in f:
            genes.append(line.strip())
    return genes
