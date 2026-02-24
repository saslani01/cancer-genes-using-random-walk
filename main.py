from loader import *
from analyzer import *
import numpy as np
import os
import time

def main():
    os.makedirs("outputs", exist_ok=True)
    start = time.time()

    # fill the adj matrix
    adj_matrix, protein_index = protein_interaction_adj_matrix("./data/interacting_proteins.txt", True)
    
    # onco genes
    onco_genes = cancer_genes("./data/onco_genes.txt")
    onco_genes_indices = [protein_index[onco_gene] for onco_gene in onco_genes]
    
    for i in range(2):
        label = "Randomized" if i == 1 else "Normal"
        # randomize 20 times
        adj_matrix = edge_swap_randomize(adj_matrix, 20 * np.count_nonzero(np.triu(adj_matrix, k=1))) if i == 1 else adj_matrix 

        print(f"Experiment on {label} Graph\n{"-" * (len(label) + 20)}\n")
        
        # shortest path
        print("* shorted_path_metric")
        observed_shortest_path_avg = avg_shortest_path(adj_matrix, onco_genes_indices) 
        shortest_path_means = avg_shortest_path_background_distribution(adj_matrix, protein_index, 1000, len(onco_genes_indices), onco_genes_indices)
        significance, p_val = avg_shortest_path_significance(observed_shortest_path_avg, shortest_path_means)
        print(f"significance={significance}\np_val={p_val}\n")
        if i == 0:
            background_distribution_plot(shortest_path_means, observed_shortest_path_avg, f"Dijkstra - {label}", f"outputs/dijkstra_background_{label}.png")

        # rwr
        print("* rwr_metric")
        stationary_probs = rwr(adj_matrix, onco_genes_indices, 0.3, 1000, 1e-10)
        observed_rwr_avg = avg_rwr(stationary_probs, onco_genes_indices)
        rwr_avgs = rwr_background_distribution(adj_matrix, protein_index, 1000, len(onco_genes), onco_genes_indices, 0.3, 1000, 1e-10)
        significance, p_val = rwr_significance(observed_rwr_avg, rwr_avgs)
        print(f"significance={significance}\np_val={p_val}\n")
        if i == 0:
            background_distribution_plot(rwr_avgs, observed_rwr_avg, f"RWR - {label}", f"outputs/rwr_background_{label}.png")

        if i == 0:
            rwr_mostly_visited_proteins(adj_matrix, stationary_probs, protein_index, 30 + len(onco_genes), True, True, "outputs/Sahand_Aslani_onco_predictions.png", "outputs/Sahand_Aslani_onco_predictions.txt")


    end = time.time()
    print(f"EXECUTION_TIME={end - start} s")


if __name__ == "__main__":
    main()
