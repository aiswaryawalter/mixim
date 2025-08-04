
from sklearn.cluster import KMeans
import numpy as np

# LARMix Routing (τ + Balancing)
def sample_latency_aware_path(simulation, latency_matrix, tau):
    """
    Selects a path using LARMix with bias parameter tau.
    latency_matrix: dict {(node_i, node_j): latency_value}
    tau: 0 (deterministic) -> 1 (uniform)
    """
    path = []
    current_layer = 0
    current_node = random.choice(simulation.layers[0])  # Entry node
    path.append(current_node)

    for layer in range(1, simulation.n_layers):
        next_layer_nodes = simulation.layers[layer]
        latencies = [latency_matrix[(current_node, node)] for node in next_layer_nodes]

        # Rank nodes by latency (0 = lowest latency)
        sorted_indices = np.argsort(latencies)
        rank_map = {node: rank for rank, node in enumerate([next_layer_nodes[i] for i in sorted_indices])}

        # Compute weights using Eq (1) from paper
        weights = []
        for node in next_layer_nodes:
            Ri = rank_map[node]
            lij = latency_matrix[(current_node, node)]
            weight = ((1/np.e)**(Ri*(1-tau))) * ((1/lij)**(1-tau))
            weights.append(weight)
        weights = np.array(weights) / np.sum(weights)

        next_node = np.random.choice(next_layer_nodes, p=weights)
        path.append(next_node)
        current_node = next_node

    return path
# Balancing Algorithm (Greedy Example)
def greedy_balance(scattering_matrix):
    """
    scattering_matrix: W x W numpy array for one layer
    Returns balanced routing weights preserving latency bias as much as possible.
    """
    tolerance = 1e-3
    balanced = False
    while not balanced:
        load = np.sum(scattering_matrix, axis=0)  # traffic per next-layer node
        overloaded = np.where(load > 1 + tolerance)[0]
        underloaded = np.where(load < 1 - tolerance)[0]

        if len(overloaded) == 0 and len(underloaded) == 0:
            balanced = True
            break

        # Scale down overloaded
        for idx in overloaded:
            scale = 1 / load[idx]
            scattering_matrix[:, idx] *= scale

        # Redistribute leftover to underloaded proportionally
        load = np.sum(scattering_matrix, axis=0)
        deficit = 1 - load[underloaded]
        total_deficit = np.sum(deficit)
        if total_deficit > 0:
            add_fraction = deficit / total_deficit
            scattering_matrix[:, underloaded] += (add_fraction * np.sum(load[overloaded] - 1))

    return scattering_matrix


# K-Means Layer Diversification
def diversify_layers(nodes, coords, n_layers):
    """
    nodes: list of node IDs
    coords: dict {node_id: (lat, lon)}
    n_layers: number of layers
    """
    X = np.array([coords[node] for node in nodes])
    kmeans = KMeans(n_clusters=n_layers).fit(X)
    labels = kmeans.labels_

    # Ensure diversity: each layer picks one from each cluster
    layers = [[] for _ in range(n_layers)]
    for cluster_id in range(n_layers):
        cluster_nodes = [node for node, label in zip(nodes, labels) if label == cluster_id]
        for i, node in enumerate(cluster_nodes):
            layer_index = i % n_layers
            layers[layer_index].append(node)

    return layers
