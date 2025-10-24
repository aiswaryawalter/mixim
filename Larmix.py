import numpy as np
import random
from sklearn.cluster import KMeans


class Larmix:
    """
    LARMix (Latency-Aware Routing for Mix networks) implementation
    Handles all latency-aware path selection and network diversification
    """
    
    def __init__(self, simulation, network):
        self.simulation = simulation
        self.network = network
        self.latency_matrix = simulation.latency_matrix
        self.node_coords = simulation.node_coords
    
    def sample_larmix_path_stratified(self, tau, sender_server_id, receiver_server_id):
        """
        Latency-aware path selection for stratified topology.
        One node per layer chosen based on latency and tau parameter.
        """
        path = []
        network_dict = self.network.network_dict
        num_layers = self.simulation.n_layers

        # Find entry mix with lowest latency from sender
        entry_candidates = network_dict[1]
        entry_latencies = []

        for node in entry_candidates:
            key = (sender_server_id, node.server_id)
            if key in self.latency_matrix:
                latency = self.latency_matrix[key]
            else:
                latency = 50.0  # default
                print(f"[DEBUG] Missing latency for {key}, using default: {latency}")
            entry_latencies.append(latency)

        # Layer 1 Weight by inverse latency (lower latency = higher weight)
        entry_weights = 1 / (np.array(entry_latencies) + 1e-6)
        entry_weights = entry_weights / np.sum(entry_weights)
        current_node = np.random.choice(entry_candidates, p=entry_weights)
        path.append(current_node)
        
        print(f"[LARMix] Layer 1 -- nodeID: {current_node.id} -- serverID: {current_node.server_id}")

        # Pick one node from each remaining layer
        for layer_idx in range(2, num_layers):
            next_layer_nodes = network_dict[layer_idx]

            # Get latencies
            latencies = [
                self.latency_matrix.get((current_node.server_id, node.server_id), 50.0)
                for node in next_layer_nodes
            ]

            # Rank nodes based on latency
            sorted_indices = np.argsort(latencies)
            rank_map = {next_layer_nodes[idx]: rank for rank, idx in enumerate(sorted_indices)}
            
            # Compute weights
            weights = []
            for node in next_layer_nodes:
                rank = rank_map[node]
                lij = self.latency_matrix.get((current_node.server_id, node.server_id), 50.0)
                weight = ((1 / np.e) ** (rank * (1 - tau))) * ((1 / lij) ** (1 - tau))
                weights.append(weight)

            weights = np.array(weights) / np.sum(weights)
            next_node = np.random.choice(next_layer_nodes, p=weights)

            path.append(next_node)
            current_node = next_node
            print(f"[LARMix] Layer: {layer_idx} -- nodeID: {next_node.id} -- serverID: {next_node.server_id}")
        
        # Select exit mix considering both current node latency AND receiver location
        if num_layers > 1:
            exit_candidates = network_dict[num_layers]
            exit_latencies = []
            
            for node in exit_candidates:
                # Combined latency: current -> exit + exit -> receiver
                current_to_exit = self.latency_matrix.get((current_node.server_id, node.server_id), 0.05)
                exit_to_receiver = self.latency_matrix.get((node.server_id, receiver_server_id), 0.05)
                total_latency = current_to_exit + exit_to_receiver
                exit_latencies.append(total_latency)
            
            # Apply tau-based weighting for exit selection
            sorted_indices = np.argsort(exit_latencies)
            rank_map = {exit_candidates[idx]: rank for rank, idx in enumerate(sorted_indices)}
            
            exit_weights = []
            for node in exit_candidates:
                rank = rank_map[node]
                total_latency = exit_latencies[exit_candidates.index(node)]
                weight = ((1 / np.e) ** (rank * (1 - tau))) * ((1 / total_latency) ** (1 - tau))
                exit_weights.append(weight)
            
            exit_weights = np.array(exit_weights) / np.sum(exit_weights)
            exit_node = np.random.choice(exit_candidates, p=exit_weights)
            path.append(exit_node)
            
            print(f"[LARMix] Last Layer -- nodeID: {exit_node.id} -- serverID: {exit_node.server_id})")

        print(f"[LARMix] Latency Aware Path: {path}")
        return path

    def sample_larmix_path_free_route(self, tau, sender_server_id, receiver_server_id, n_hops):
        """
        Sample a latency-aware path through free route network
        """
        # Get all mixes (they're all in layer 1 for free route)
        all_mixes = self.network.network_dict[1] if 1 in self.network.network_dict else []
        
        if len(all_mixes) < n_hops:
            print(f"[WARNING] Not enough mixes for {n_hops} hops, using {len(all_mixes)} mixes")
            n_hops = len(all_mixes)
        
        # Start from sender location
        current_server_id = sender_server_id
        path = []
        used_mixes = set()
        
        print(f"[LARMix Free Route] Starting path from server {current_server_id}")
        
        for hop in range(n_hops):
            # Calculate latencies from current position to all unused mixes
            candidate_mixes = [mix for mix in all_mixes if mix not in used_mixes]
            
            if not candidate_mixes:
                break
                
            latencies = []
            for mix in candidate_mixes:
                mix_server_id = int(mix.server_id)
                latency = self.latency_matrix.get((current_server_id, mix_server_id), 0.1)
                latencies.append(latency)
            
            # Apply tau-based selection (Boltzmann exploration)
            if tau == 0:
                # Deterministic: choose lowest latency
                best_idx = np.argmin(latencies)
                chosen_mix = candidate_mixes[best_idx]
            elif tau >= 1:
                # Uniform random
                chosen_mix = np.random.choice(candidate_mixes)
            else:
                # Boltzmann exploration
                inv_latencies = [-lat/tau for lat in latencies]  # Negative because we want low latency
                exp_values = np.exp(inv_latencies)
                probabilities = exp_values / np.sum(exp_values)
                chosen_mix = np.random.choice(candidate_mixes, p=probabilities)
            
            path.append(chosen_mix)
            used_mixes.add(chosen_mix)
            current_server_id = int(chosen_mix.server_id)
            
            print(f"[LARMix Free Route] Hop {hop+1}: Mix {chosen_mix.id} (Server {current_server_id})")
        
        # Calculate final latency to receiver
        final_latency = self.latency_matrix.get((current_server_id, receiver_server_id), 0.1)
        print(f"[LARMix Free Route] Final hop to receiver: latency {final_latency}")
        
        return path

    def diversify_layers(self, mix_ids, coords, num_layers, mixes_per_layer):
        """
        Creates geographically diversified layers for stratified topology.
        """
        # Prepare coordinate matrix for clustering
        X = np.array([coords[mix_id] for mix_id in mix_ids])
        kmeans = KMeans(n_clusters=mixes_per_layer, random_state=42).fit(X)
        labels = kmeans.labels_

        # Group mixes by cluster
        cluster_to_nodes = {i: [] for i in range(mixes_per_layer)}
        for mix_id, label in zip(mix_ids, labels):
            cluster_to_nodes[label].append(mix_id)

        # Create layers by picking one node from each cluster for every layer
        layers = [[] for _ in range(num_layers)]
        used_nodes = set()
        for layer_idx in range(num_layers):
            for cluster_idx in range(mixes_per_layer):
                candidates = [node for node in cluster_to_nodes[cluster_idx] if node not in used_nodes]
                if not candidates:  # cluster exhausted → pick any unused node
                    remaining = [node for node in mix_ids if node not in used_nodes]
                    if not remaining:
                        raise ValueError("Not enough unique mixes to fill layers.")
                    node = random.choice(remaining)
                else:
                    node = random.choice(candidates)
                layers[layer_idx].append(node)
                used_nodes.add(node)

        return layers

    def greedy_balance(self, scattering_matrix):
        """
        Greedy load balancing algorithm for LARMix
        """
        tolerance = 1e-3
        balanced = False
        while not balanced:
            load = np.sum(scattering_matrix, axis=0)
            overloaded = np.where(load > 1 + tolerance)[0]
            underloaded = np.where(load < 1 - tolerance)[0]
            if len(overloaded) == 0 and len(underloaded) == 0:
                balanced = True
                break
            for idx in overloaded:
                scale = 1 / load[idx]
                scattering_matrix[:, idx] *= scale
            load = np.sum(scattering_matrix, axis=0)
            deficit = 1 - load[underloaded]
            total_deficit = np.sum(deficit)
            if total_deficit > 0:
                add_fraction = deficit / total_deficit
                scattering_matrix[:, underloaded] += (add_fraction * np.sum(load[overloaded] - 1))
        return scattering_matrix

    def assign_server_ids_to_mixes(self, server_info, mixes):
        """
        Assign server IDs to mixes using round-robin approach
        """
        for i, mix in enumerate(mixes):
            server_row = server_info.iloc[i % len(server_info)]
            mix.server_id = server_row['id']
            print(f"[LARMix] Mix {mix.id} assigned to server {mix.server_id}")