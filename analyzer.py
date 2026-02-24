from scipy.sparse.csgraph import dijkstra
import numpy as np
import random 
import networkx as nx
import matplotlib.pyplot as plt

# given a selected list of protein indices, for each pair of indices calculate the shortest path using dijkstra and then return the mean
def avg_shortest_path(adj_matrix, target_protein_indices):
    distance_matrix = dijkstra(adj_matrix, indices=target_protein_indices)
    
    distances = []
    for i in range(len(target_protein_indices)):
        for j in range(i + 1, len(target_protein_indices)):
            distances.append(distance_matrix[i][target_protein_indices[j]])

    return np.mean(distances)

def avg_shortest_path_background_distribution(adj_matrix, protein_index, trial_count, n, exclude_indices):
    mean_distances = []
    non_cancer_protein_indices = [i for i in protein_index.values() if i not in exclude_indices]       
    for i in range(trial_count):
        # random non-cancer proteins sampling 
        target_protein_indices = random.sample(non_cancer_protein_indices, n)
        mean_distance = avg_shortest_path(adj_matrix, target_protein_indices)
        mean_distances.append(mean_distance)
        #print(f"trial-{i}={mean_distance}")
    return mean_distances

def avg_shortest_path_significance(observed_shortest_path_mean, shortest_path_means, alpha=0.1):
    count = sum(1 for shortest_path_mean in shortest_path_means if shortest_path_mean <= observed_shortest_path_mean) 
    p_val = count / len(shortest_path_means)
    
    return p_val < alpha, p_val

def rwr(adj_matrix, start_protein_indices, gamma, max_iterations, thresh):
    protein_count = adj_matrix.shape[0]
    
    neighbor_count = adj_matrix.sum(axis=1)
    
    transition_matrix = adj_matrix / neighbor_count.reshape(-1, 1)
    
    start_prob = np.zeros(protein_count)
    for p_i in start_protein_indices:
        start_prob[p_i] = 1 / len(start_protein_indices)
    curr_prob = start_prob.copy()
    
    for i in range(max_iterations):
        walk_prob = transition_matrix.T @ curr_prob
        next_prob = (1 - gamma) * walk_prob + gamma * start_prob

        change = np.linalg.norm(next_prob - curr_prob)
        curr_prob = next_prob

        #print(f"magnitude-difference={change}")
        if change < thresh:
            #print(f"[CONVERGED] i={i + 1} \n")
            break
    else:
        print(f"[NO CONVERGENCE] (max_iteration={max_iterations})")
        exit(1)


    return curr_prob

def avg_rwr(stationary_probs, onco_gene_indices):
    return np.mean(stationary_probs[onco_gene_indices])

def rwr_background_distribution(adj_matrix, protein_index, trial_count, n, exclude_indices, gamma, max_iterations, thresh):
    non_cancer_protein_indices = [i for i in protein_index.values() if i not in exclude_indices]
    rwr_avgs = []
    for i in range(trial_count):
        sample_indices = random.sample(non_cancer_protein_indices, n)
        stationary_probs = rwr(adj_matrix, sample_indices, gamma, max_iterations, thresh)
        rwr_avg = avg_rwr(stationary_probs, sample_indices)
        rwr_avgs.append(rwr_avg)
        #print(f"trial-{i}={rwr_avg}")
    return rwr_avgs

def rwr_significance(observed_rwr_avg, rwr_avgs, alpha=0.10):
    count = sum(1 for rwr_avg in rwr_avgs if rwr_avg >= observed_rwr_avg)
    p_val = count / len(rwr_avgs)
    
    return p_val < alpha, p_val

# k mostly visited proteins
def rwr_mostly_visited_proteins(adj_matrix, stationary_probs, protein_index, k, vizualize, save_top_k, visual_path, txt_path):
    protein_index_flipped  = {i: name for name, i in protein_index.items()}
        
    # top k
    indices = np.argsort(stationary_probs)[::-1][:k]
    names = [protein_index_flipped[i] for i in indices]
    stationary_probs = stationary_probs[indices] # no mutaion; we are good

    if vizualize:
        graph = nx.from_numpy_array(adj_matrix[np.ix_(indices, indices)])
        graph = nx.relabel_nodes(graph, {i: f"{name}\n{prob * 100:.2f}" for i, (name, prob) in enumerate(zip(names, stationary_probs))})        
        colors = [p * 100 for p in stationary_probs]
        
        _, ax = plt.subplots(figsize=(15, 15))
        pos = nx.spring_layout(graph, k=2, scale=3, seed=13)
        nx.draw(graph, ax=ax, pos=pos, node_color=colors, cmap=plt.cm.Blues, node_size=500, with_labels=True)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.Blues, norm=plt.Normalize(min(colors), max(colors)))
        ax.set_title(f"Stationary Probability of {k} Mostly Visited Proteins", fontsize=20)
        cbar = plt.colorbar(sm, ax=ax)
        cbar.ax.tick_params(labelsize=13)
        cbar.set_label("Stationary Probability (%)", fontsize=13)        
        plt.savefig(visual_path)
        plt.close()

    if save_top_k:
        with open(txt_path, 'w') as f:
            for name in names:
                f.write(f"{name}\n")

    return names, stationary_probs

# edge swappigng degree preserving randomization
def edge_swap_randomize(adj_matrix, swap_count):
    graph = nx.from_numpy_array(adj_matrix)
    nx.double_edge_swap(graph, nswap=swap_count, max_tries=swap_count * 2)
    return nx.to_numpy_array(graph)

def background_distribution_plot(background_scores, observed_score, metric_name, save_path):
    _, ax = plt.subplots(figsize=(10, 6))
    ax.hist(background_scores, bins=30, color="lightblue", edgecolor="white", label="Random Samples")
    ax.axvline(observed_score, color="red", linewidth=2, label=f"Oncogene Score = {observed_score:.4f}")
    ax.set_xlabel("Average Proximity")
    ax.set_ylabel("Count")
    ax.set_title(f"Background Distribution - {metric_name}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
