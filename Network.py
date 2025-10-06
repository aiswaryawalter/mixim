from PoissonMix import PoissonMix
from Pool import Pool
from TimedMix import TimedMix

import random
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans


class Network:
    all_mixes = []
    network_dict = {}  # 1:[list of mixes in layer 1], 2:[list of mixes in layer 2], ...

    def __init__(self, mix_type, num_layers, nbr_mixes_layers, corrupt, unifrom_corruption, simulation,
                 threshold,
                 flush_percent, topology,fully_connected, flushtime, probability_dist_mixes, n_cascades, 
                 m_barabasi_mixes,  link_based_dummies, multiple_hop_dummies, rate_mix_dummies, Network_template, numberTargets):
        self.simulation = simulation
        self.num_layers = num_layers
        self.mix_type = mix_type
        self.mixesPerLayer = nbr_mixes_layers
        self.corrupt = corrupt
        self.unifrom_corruption = unifrom_corruption
        self.env = simulation.env
        self.threshold = threshold
        self.flush_percent = flush_percent
        self.topology = topology
        self.fully_connected = fully_connected
        self.flushtime = flushtime
        self.n_cascades = n_cascades
        self.m_barabasi_mixes = m_barabasi_mixes
        self.link_based_dummies = link_based_dummies
        self.multiple_hop_dummies =multiple_hop_dummies
        self.rate_mix_dummies = rate_mix_dummies
        self.probability_dist_mixes = probability_dist_mixes
        self.Network_template = Network_template
        self.numberTargets = numberTargets
        self.list_cascades = {}
        self.n_cascades = 6

        # larmix
        self.latency_matrix = simulation.latency_matrix
        self.node_coords = simulation.node_coords

        self.server_info = pd.read_csv('servers.csv')

        self.create_network()

    def create_network(self):
        mixnb = 1
        self.all_mixes = set()
        if self.topology == 'stratified':
            Nbr_Corruption = 0
            for layer in range(1, self.num_layers + 1):
                c = 0
                self.network_dict[layer] = []
                for _ in range(self.mixesPerLayer):
                    if self.unifrom_corruption:
                        if c < self.corrupt / self.simulation.n_layers:
                            varCorrupt = True
                            c += 1
                        else:
                            varCorrupt = False
                    else:
                        if Nbr_Corruption < self.corrupt:
                            varCorrupt = random.choice([True, False])
                            if varCorrupt:
                                Nbr_Corruption += 1
                        else:
                            varCorrupt = False
                    mix = self.get_mixnode(self.mix_type, mixnb, layer, self.numberTargets, varCorrupt,
                                           self.probability_dist_mixes[layer - 1][_])
                    # larmix
                    server_row = self.server_info.iloc[(mixnb - 1) % len(self.server_info)]
                    mix.server_id = server_row['id']
                    print(f"[DEBUG] Created mix - Node ID: {mix.id}, Server ID: {mix.server_id}, Layer: {mix.layer}")

                    # mix.location = server_row['location']
                    self.all_mixes.add(mix)
                    self.network_dict[layer] += [mix]
                    mixnb += 1
            if self.simulation.routing == "larmix":
                # prepare coordinates mapping
                coords = {row['id']: (row['latitude'], row['longitude'])
                        for _, row in self.server_info.iterrows()}

                mix_ids = [mix.server_id for mix in self.all_mixes]
                print(f"[Larmix] Mix ids ==> {mix_ids}")
                print(f"[DEBUG] Before diversification:")
                for mix in self.all_mixes:
                    print(f"  Mix ID: {mix.id}, Server ID: {mix.server_id}, Layer: {mix.layer}")
                diversified_layers = self.diversify_layers(mix_ids, coords, self.num_layers, self.mixesPerLayer)
                if any(len(layer) == 0 for layer in diversified_layers):
                    print("[ERROR] Diversification produced empty layer.")
                    return

                # rebuild network_dict based on diversified layers
                self.network_dict = {}
                for layer_id, layer_server_ids in enumerate(diversified_layers, start=1):
                    self.network_dict[layer_id] = [mix for mix in self.all_mixes
                                                    if mix.server_id in layer_server_ids]
                    for mix in self.network_dict[layer_id]:
                        mix.layer = layer_id
                    
                    print(f"[DEBUG] Layer {layer_id} after diversification:")
                    for mix in self.network_dict[layer_id]:
                        print(f"  Mix ID: {mix.id}, Server ID: {mix.server_id}, New Layer: {mix.layer}")
                    # neighbor assignment
                    if layer_id < self.num_layers and (layer_id + 1) in self.network_dict:
                        for mix in self.network_dict[layer_id]:
                            mix.neighbors = self.network_dict[layer_id + 1]
                
            else:
                for mix in self.all_mixes:
                    if self.fully_connected:
                        if mix.layer + 1 in self.network_dict:  # last mix doesn't need neighbors
                            mix.neighbors = self.network_dict[mix.layer + 1]
                        if mix.layer == self.simulation.n_layers:
                            mix.neighbors = self.network_dict[1]
                    else:
                        pass
                        #for mix in self.MixesAll:
                            #if mix.id == 1:
                                #mix.neighbors = []
                                #mix.neighbors.append(self.LayerDict[mix.layer + 1][0])
                                #mix.neighbors.append(self.LayerDict[mix.layer + 1][1])
        elif self.topology == 'cyclic_stratified':
            for layer in range(1, self.simulation.n_layers + 1):
                self.network_dict[layer] = []
                for _ in range(self.simulation.n_mixes_per_layer):
                    varCorrupt = False
                    mix = self.get_mixnode(
                        self.mix_type,
                        mixnb,
                        layer,
                        self.numberTargets,
                        varCorrupt,
                        self.probability_dist_mixes[layer - 1][_]
                    )
                    self.all_mixes.add(mix)
                    self.network_dict[layer].append(mix)
                    mixnb += 1
                    # larmix
                    # server_row = self.server_info.iloc[(mixnb - 1) % len(self.server_info)]
                    # mix.server_id = server_row['server_id']
            for layer in range(1, self.simulation.n_layers + 1):
                next_layer = (layer % self.simulation.n_layers) + 1
                for mix in self.network_dict[layer]:
                    mix.neighbors = self.network_dict[next_layer]

        elif self.topology == 'XRD':
            mixnb = 1
            for n in range(1, 1 + self.n_cascades):
                cascade = []
                for m in range(self.num_layers):
                    varCorrupt = False
                    mix = self.get_mixnode(self.mix_type, mixnb, m + 1, self.numberTargets, varCorrupt,
                                           1 / self.n_cascades)
                    mix.n_chain = n
                    self.all_mixes.add(mix)
                    mixnb += 1
                    cascade.append(mix)
                self.list_cascades[n] = cascade
            for n, list in self.list_cascades.items():
                print('Chain number', n, ':', list)
        
        elif self.topology == 'free route':
            self.network_dict[1] = []  # if we treat everything as "layer 1"
            self.all_mixes = set()
            Nbr_Corruption = 0
            
            for i in range(self.mixesPerLayer):  
                if self.unifrom_corruption:
                    # (Same logic as 'stratified' to spread corruption)
                    varCorrupt = (Nbr_Corruption < self.corrupt)
                    if varCorrupt:
                        Nbr_Corruption += 1
                else:
                    # or pick randomly until you reach self.corrupt
                    if Nbr_Corruption < self.corrupt:
                        varCorrupt = random.choice([True, False])
                        if varCorrupt:
                            Nbr_Corruption += 1
                    else:
                        varCorrupt = False
                # Create the mix (poisson, timed, or pool) exactly like stratified:
                mix = self.get_mixnode(
                    self.mix_type,
                    i+1,                 # mix ID
                    1,                   # position = 1 if no layers
                    self.numberTargets,
                    varCorrupt,
                    weight_mix=1.0/self.mixesPerLayer
                )
                self.all_mixes.add(mix)
                self.network_dict[1].append(mix)
                # larmix
                # server_row = self.server_info.iloc[i % len(self.server_info)]
                # mix.server_id = server_row['server_id']


            # Define a connectivity or neighbor relationship for "free route"
            list_of_mixes = self.network_dict[1]
            for mix in list_of_mixes:
                # Suppose we want each mix to have k random neighbors
                # or all the other mixes if "fully_connected=True"
                if self.fully_connected:
                    mix.neighbors = [m for m in list_of_mixes if m != mix]
                else:
                    # for partial connectivity, pick some random subset:
                    possible_neighbors = [m for m in list_of_mixes if m != mix]
                    # pick e.g. 2 random neighbors:
                    mix.neighbors = random.sample(possible_neighbors, k=2)
       
        elif self.topology == 'ba topology':
            N = self.mixesPerLayer  
            m = self.m_barabasi_mixes 
            # build adjacency list
            adjacency_list = self.ba_adjacency(N, m)

            # create mixes
            self.network_dict[1] = []
            for node_id in range(N):
                varCorrupt = False
                mix = self.get_mixnode(
                    self.mix_type,
                    id=node_id + 1,        
                    position=1,          
                    numberTargets=self.numberTargets,
                    corrupt=varCorrupt,
                    weight_mix=1.0 / N     
                )
                self.network_dict[1].append(mix)
                self.all_mixes.add(mix)
                # larmix
                # server_row = self.server_info.iloc[i % len(self.server_info)]
                # mix.server_id = server_row['server_id']
                
            # assign neighbors based on adjacency_list
            for node_id, mix in enumerate(self.network_dict[1]):
                neighbor_ids = adjacency_list[node_id]
                # Convert neighbor indices to actual mix objects
                mix.neighbors = [self.network_dict[1][nbr_id] for nbr_id in neighbor_ids]



    def get_mixnode(self, mix_type, id, position, numberTargets, corrupt, weight_mix):
        if mix_type == 'poisson':
            return PoissonMix(id, self.simulation, position, self.link_based_dummies,self.multiple_hop_dummies, self.rate_mix_dummies,
                              numberTargets, corrupt, weight_mix)
        elif mix_type == 'pool':
            return Pool(id, self.simulation, position, self.threshold, self.flush_percent, numberTargets,
                        corrupt, weight_mix)
        elif mix_type == 'time':
            return TimedMix(id, self.simulation, position, self.flushtime, numberTargets, corrupt,
                            weight_mix)

    def odd(self, number):
        return number % 2 == 1
    
    def ba_adjacency(self, N, m):
        if m < 1 or m >= N:
            raise ValueError("BA-topology parameter m must be in [1, N-1].")

        # adjacency_list[i] = list of neighbors of node i
        adjacency_list = [[] for _ in range(N)]

        # fully connected (complete) graph of m nodes
        for i in range(m):
            for j in range(i+1, m):
                adjacency_list[i].append(j)
                adjacency_list[j].append(i)

        # degree array for each node
        degree = [0]*N
        for i in range(m):
            degree[i] = m-1 

        # add remaining nodes one by one
        for new_node in range(m, N):
            degree_sum = sum(degree[:new_node])

            # pick m distinct existing nodes using preferential attachment
            connected = set()
            while len(connected) < m:
                candidate = random.randrange(new_node)
                p_attach = degree[candidate] / degree_sum
                if random.random() < p_attach:
                    connected.add(candidate)

            # connect new_node with each in connected
            for cand in connected:
                adjacency_list[new_node].append(cand)
                adjacency_list[cand].append(new_node)
                degree[new_node] += 1
                degree[cand] += 1

        return adjacency_list

    def greedy_balance(scattering_matrix):
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
    
    def diversify_layers(self, mix_ids, coords, num_layers, mixes_per_layer):
        """
        Creates geographically diversified layers.
        
        Args:
            mix_ids: list of mix server_ids (one per mix node)
            coords: dict {server_id: (lat, lon)}
            num_layers: number of layers (stratified depth)
            mixes_per_layer: number of mixes per layer
        
        Returns:
            layers: list of layers, each layer = list of server_ids
        """
        # --- prepare coordinate matrix for clustering ---
        X = np.array([coords[mix_id] for mix_id in mix_ids])
        kmeans = KMeans(n_clusters=mixes_per_layer, random_state=42).fit(X)
        labels = kmeans.labels_

        # --- group mixes by cluster ---
        cluster_to_nodes = {i: [] for i in range(mixes_per_layer)}
        for mix_id, label in zip(mix_ids, labels):
            cluster_to_nodes[label].append(mix_id)

        # --- create layers by picking one node from each cluster for every layer ---
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

    def sample_latency_aware_path(self, tau, sender_server_id, receiver_server_id):
        """
        Latency-aware path selection for stratified topology.
        One node per layer chosen based on latency and tau parameter.
        """
        path = []

        # Find entry mix with lowest latency from sender
        entry_candidates = self.network_dict[1]
        entry_latencies = []
        # print(f"[DEBUG] Sender server ID: {sender_server_id} (type: {type(sender_server_id)})")
        # print(f"[DEBUG] Available latency matrix keys sample: {list(self.latency_matrix.keys())[:5]}")

        for node in entry_candidates:
            key = (sender_server_id, node.server_id)
            if key in self.latency_matrix:
                latency = self.latency_matrix[key]
                # print(f"[DEBUG] Found latency for {key}: {latency}")
            else:
                latency = 50.0  # default
                print(f"[DEBUG] Missing latency for {key}, using default: {latency}")
            entry_latencies.append(latency)

        
        # Layer 1 Weight by inverse latency (lower latency = higher weight)
        entry_weights = 1 / (np.array(entry_latencies) + 1e-6)
        entry_weights = entry_weights / np.sum(entry_weights)
        current_node = np.random.choice(entry_candidates, p=entry_weights)
        path.append(current_node)
        
        # print(f"[LARMix] Layer 1 -- nodeID: {current_node.id} -- serverID: {current_node.server_id}")

        # --- Pick one node from each remaining layer ---
        for layer_idx in range(2, self.num_layers):
            next_layer_nodes = self.network_dict[layer_idx]

            # Get latencies
            latencies = [
                self.latency_matrix.get((current_node.server_id, node.server_id), 50.0)
                for node in next_layer_nodes
            ]
            # print(f"[LARMix] Latencies from ServerID {current_node.server_id}: {latencies}")

            # Rank nodes based on latency
            sorted_indices = np.argsort(latencies)
            rank_map = {next_layer_nodes[idx]: rank for rank, idx in enumerate(sorted_indices)}

            # print(f"[DEBUG] Sorted_indices of Layer {layer_idx}: {sorted_indices}")
            # print(f"[DEBUG] Rank map of Layer {layer_idx}: {rank_map}")
            
            # Compute weights
            weights = []
            for node in next_layer_nodes:
                rank = rank_map[node]
                lij = self.latency_matrix.get((current_node.server_id, node.server_id), 50.0)
                weight = ((1 / np.e) ** (rank * (1 - tau))) * ((1 / lij) ** (1 - tau))
                # print(f"[DEBUG] Mixnode {node.id} -> from Server ID {current_node.server_id} to {node.server_id} Latency : {lij}, Weight: {weight}")
                weights.append(weight)

            weights = np.array(weights) / np.sum(weights)
            # print(f"[DEBUG] Normalized weights: {weights}")
            next_node = np.random.choice(next_layer_nodes, p=weights)

            path.append(next_node)
            current_node = next_node
            # print(f"[LARMix] Layer: {layer_idx} -- nodeID: {next_node.id} -- serverID: {next_node.server_id}")
        
        # Select exit mix considering both current node latency AND receiver location
        if self.num_layers > 1:
            exit_candidates = self.network_dict[self.num_layers]
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
            
            # print(f"[LARMix] Last Layer -- nodeID: {exit_node.id} -- serverID: {exit_node.server_id})")
            # print(f"[LARMix] Exit total latencies: {list(zip([n.server_id for n in exit_candidates], exit_latencies))}")

        # print(f"[LARMix] Latency Aware Path: {path}")

        return path
