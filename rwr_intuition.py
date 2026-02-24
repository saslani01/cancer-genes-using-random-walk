from analyzer import rwr
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import os

graph = nx.Graph([
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 2),
    (1, 3),
    (2, 0),
    (2, 1),
    (2, 3),
    (2, 4),
    (3, 1),
    (3, 2),
    (4, 5),
    (4, 2),
    (4, 6),
    (5, 4),
    (6, 4)
])

adj_matrix = nx.to_numpy_array(graph)
start_node = [0]

def main():
    os.makedirs("outputs", exist_ok=True)
    
    stationary_probs = rwr(adj_matrix, start_node, 0.3, 1000, 1e-10)
    rwr_plot(graph, stationary_probs, "Random Walk with Restart (gamma=0.3)", "outputs/rwr_simulation.png")

    adj_weighted = adj_matrix.copy()
    adj_weighted[2][4] = adj_weighted[4][2] = 100
    adj_weighted[0][1] = adj_weighted[1][0] = 0.01

    stationary_probs = rwr(adj_weighted, start_node, 0.3, 1000, 1e-10)
    rwr_plot(graph, stationary_probs, "Random Walk with Restart Weighted (gamma=0.3)", "outputs/rwr_weighted_simulation.png")

    rwr_plot(graph, rwr(adj_matrix, start_node, 0.0, 1000, 1e-10), "Random Walk with Restart (gamma=0)", "outputs/rwr_simulation_gamma_0.png")
    rwr_plot(graph, rwr(adj_matrix, start_node, 1.0, 1000, 1e-10), "Random Walk with Restart (gamma=1)", "outputs/rwr_simulation_gamma_1.png")

    print("Spearman Correlation between Degrees and Stationary Probability")
    gammas = [1.0, 0.0, 0.3]
    degrees = [d for _, d in graph.degree()]
    for gamma in gammas:
        stationary_probs = rwr(adj_matrix, start_node, gamma, 1000, 1e-10)
        corr, p_val = spearmanr(degrees, stationary_probs)
        print(f"spearmam_corr={corr:.4f}, p={p_val:.4f}, gamma={gamma}")

def rwr_plot(G, stationary_probs, title, path):
    pos = nx.spring_layout(G, seed=42)
    labels = {i: f"{i}\n{p * 100:.2f}%" for i, p in enumerate(stationary_probs)}
    _, ax = plt.subplots(figsize=(8, 6))
    nx.draw(G, pos=pos, ax=ax, node_color=list(stationary_probs), cmap=plt.cm.Blues, node_size=1000, labels=labels, with_labels=True)
    ax.set_title(title, fontsize=14)
    plt.savefig(path)
    plt.close()

if __name__ == "__main__":
    main()