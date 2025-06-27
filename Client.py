from random import choice, sample
from Message import Message
from numpy.random import exponential
import numpy as np
import random

class Client:
    def __init__(self, simulation, id, network_dict, rate_client, mu, probability_dist_mixes, n_targets, n_hops, client_dummies, rate_client_dummies, Log):
        self.id = id
        self.env = simulation.env
        self.simulation = simulation  # simulation object
        self.network_dict = network_dict
        self.class_ends = self.env.event()
        self.mu = mu  # Avg delay: for delays at the poisson mixes
        self.probability_dist_mixes = probability_dist_mixes
        self.other_clients = set()
        self.message_id = 1

        self.rate_client = rate_client
        self.all_mixes = []
        self.n_targets = n_targets
        self.n_hops = n_hops
        self.client_dummies = client_dummies
        self.rate_client_dummies = rate_client_dummies
        self.log = Log
        if self.simulation.topology == 'stratified':
            for layer in range(1, len(self.network_dict) + 1):
                self.all_mixes += self.network_dict[layer]
        elif self.simulation.topology == 'ba topology':
            for layer in range(1, len(self.network_dict) + 1):
                self.all_mixes += self.network_dict[layer]
        elif self.simulation.topology == 'XRD':
            self.set_chains = self.network_dict
        if self.simulation.topology == 'free route':
            for layer in range(1, len(self.network_dict) + 1):
                self.all_mixes += self.network_dict[layer]
        self.env.process(self.send_message('Real', self.rate_client))
        if self.client_dummies:
            self.env.process(self.send_message('ClientDummy', self.rate_client_dummies))

    def create_message(self, message_type, rate_client):
        np.random.seed()
        delay_client = exponential(rate_client)
        route = [self]
        route_ids = [self.id]
        delays = [delay_client]
        pr_target = [0.0 for _ in range(self.n_targets)]
        if (self.simulation.topology == 'free route' and 
            self.simulation.routing == 'source'):
            # self.n_hops = np.random.randint(2, 6)  
            for _ in range(self.n_hops):
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
                node = choice(self.all_mixes)  
                while node in route:
                    node = choice(self.all_mixes)
                route.append(node)
                route_ids.append(node.id)
            print(f"[Free Route Debug] Route: {route}")
            
        elif (self.simulation.topology == 'ba topology' and 
            self.simulation.routing == 'source'):
        
            current_node = random.choice(self.all_mixes)

            # Append the first node to the route
            route.append(current_node)
            route_ids.append(current_node.id)

            # For each hop in the path, choose one random neighbor
            for _ in range(self.n_hops - 1):
                neighbors = list(current_node.neighbors)
                if not neighbors:
                    break

                node_next = random.choice(neighbors)
                
                # avoid repeating nodes
                while node_next in route:
                    node_next = random.choice(neighbors)
                
                # Append a delay for this hop
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
                
                # Append the next node to the route
                route.append(node_next)
                route_ids.append(node_next.id)

                current_node = node_next

            print(f"[BA Debug] route so far: {route}")  

        elif (self.simulation.topology == 'cyclic_stratified' and 
            self.simulation.routing == 'source'):
            start_layer = random.randint(1, self.simulation.n_layers)
            print(f"[Debug] Start Layer: {start_layer}")  
            current_layer = start_layer
            prev_node     = None 

            for layer in range(self.simulation.n_layers):
                # choose the mix for the *current* layer
                if prev_node is None:                    # first hop
                    node = np.random.choice(
                        self.network_dict[current_layer],
                        p=self.probability_dist_mixes[current_layer - 1]  # correct index
                    )
                else:                                   # subsequent hops
                    node = np.random.choice(prev_node.neighbors)

                # record hop‑delay and route entry
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
                route.append(node)
                route_ids.append(node.id)

                # prepare for next iteration
                prev_node     = node
                current_layer = (current_layer % self.simulation.n_layers) + 1
        elif (
            self.simulation.routing == 'larmix'
            and self.simulation.topology in ['stratified']):            
            best_route = self.sample_latency_aware_path()
            route = [self] + best_route
            route_ids = [node.id for node in route]
            delays = [exponential(self.mu) for _ in range(len(route) - 1)]    
        else:
            for layer in range(1, self.simulation.n_layers+1):
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
                if self.simulation.routing == 'source' and self.simulation.topology == 'stratified'\
                        or (self.simulation.routing == 'hopbyhop' and self.simulation.topology == 'stratified' and layer == 1):
                    if self.simulation.n_layers ==1 and self.simulation.n_mixes_per_layer == 1:
                        node = self.network_dict[1][0]
                        route.append(node)
                        route_ids.append(node.id)
                    else:
                        if layer == 1:
                            node = np.random.choice(self.network_dict[layer], p=self.probability_dist_mixes[layer - 1])
                            route.append(node)
                            route_ids.append(node.id)
                        else:
                            node = np.random.choice(node.neighbors)
                            route.append(node)
                            route_ids.append(node.id)
                # elif self.simulation.routing == 'source' and self.simulation.topology == 'free route'\
                elif (self.simulation.routing == 'hopbyhop' and self.simulation.topology == 'free route' and layer == 1):
                    node = choice(self.all_mixes)
                    while node in route:
                        node = choice(self.all_mixes)
                    route.append(node)
                    route_ids.append(node.id)
                    print(f"==>> route: {route}")
                elif self.simulation.routing == 'hopbyhop' and layer != 1:
                    route.append(None)
                    route_ids.append(None)
                elif self.simulation.routing == 'source' and self.simulation.topology == 'XRD':
                    chain = random.choice(self.set_chains)
                    route = [self]
                    route_ids = [self.id]
                    for node in chain:
                        route.append(node)
                        route_ids.append((node.id))
        delays += [0]
        receiver = sample(list(self.other_clients), k=1)[0]
        print(f"[Route Delays]: {delays}")
        print(f"==>> receiver: {receiver} at time {self.env.now}")
        route += [receiver]
        route_ids += [receiver.id]
        print(f"[Debug] ==>> route: {route}") 

        message = Message(self.message_id, message_type, self, route, delays, pr_target,False)
        if self.message_id == 1 and self.id ==1:
            for i in range(len(self.probability_dist_mixes)):
                if self.simulation.printing:
                    print("Weights Layer %d %s"%(i, self.probability_dist_mixes))
                else:
                    pass
        self.message_id += 1
        return message, delay_client



    def receive_message(self, message):
        message.timeReceived = self.env.now
        self.log.received_messages_f(message)
        if message.target_bool and self.simulation.printing:
            print(f'Target message arrived at destination Client at time {self.env.now}')
        if message.type == 'Real' or message.type == 'ClientDummy':
            message.route[0].receive_ack(message)


    def send_message(self, message_type, rate_client):
        while True:
            message, sending_time = self.create_message(message_type, rate_client)
            yield self.env.timeout(sending_time)
            print(f"==>> Sending Time: {sending_time}")
            message.time_left = self.env.now
            self.log.sent_messages_f(message)
            self.env.process(self.simulation.attacker.relay(message, message.route[1]))
            print(f"==>> message.route[1]: {message.route}")

    def receive_ack(self, message):  # Message received
        pass

    def __str__(self):
        return 'Client id: {}'.format(self.id)

    def __repr__(self):
        return self.__str__()
    
    def sample_latency_aware_path(self):
        valid_paths = []
        if self.simulation.topology in ['stratified', 'cyclic_stratified']:
            num_hops = self.simulation.n_layers
        else:
            num_hops = self.simulation.n_hops

        for _ in range(self.simulation.num_path_samples):
            path = []
            current_layer = 1
            valid = True
            total_latency = 0

            for hop in range(num_hops):
                if self.simulation.topology in ['stratified', 'cyclic_stratified']:
                    candidates = self.network_dict[current_layer]
                else:
                    candidates = self.all_mixes
                node = np.random.choice(candidates)
                path.append(node)

                # Estimate latency to next node if available
                if len(path) > 1:
                    prev = path[-2]
                    if node.id in prev.latency_to_neighbors:
                        total_latency += prev.latency_to_neighbors[node.id]
                    else:
                        total_latency += 2 # assume default if unknown

                if total_latency > self.simulation.latency_bound:
                    valid = False
                    break

                if self.simulation.topology in ['stratified', 'cyclic_stratified', 'XRD']:
                    current_layer = (current_layer % self.simulation.n_layers) + 1

            if valid:
                valid_paths.append((path, total_latency))
                
        if not valid_paths:
            print(f"[Warning] No valid paths under latency bound {self.simulation.latency_bound}. Using random path.")
            # Fallback to random path
            path = []
            current_layer = 1
            for _ in range(num_hops):
                if self.simulation.topology in ['stratified', 'cyclic_stratified']:
                    node = np.random.choice(self.network_dict[current_layer])
                    current_layer = (current_layer % self.simulation.n_layers) + 1
                else:
                    node = np.random.choice(self.all_mixes)
                path.append(node)
            return path

        
        entropies = [self.route_entropy(p[0]) for p in valid_paths]
        probs = self.softmax(entropies)
        chosen_index = np.random.choice(len(valid_paths), p=probs)
        best_route = valid_paths[chosen_index][0]
        print(f"[Larmix Candidate Routes]==> {[p[0] for p in valid_paths]}")
        print(f"[Larmix Candidate Entropies]==> {entropies}")
        print(f"[Larmix Best Route]==> {best_route}")
        return best_route
    
    def route_entropy(self, route):
        entropy = 0.0
        epsilon = 1e-9  # to avoid log(0)

        if self.simulation.topology in ['stratified', 'cyclic_stratified']:
            for i, mix in enumerate(route):
                layer = (i % self.simulation.n_layers)
                mix_list = self.network_dict[layer + 1]
                p = self.probability_dist_mixes[layer][mix_list.index(mix)]
                if p > 0:
                    # entropy -= p * np.log2(p)
                    entropy += -np.log2(p + epsilon) + np.random.uniform(0, 0.01)
        else:
            # Free route / BA topology — use flat list
            mix_list = self.all_mixes
            probs = self.probability_dist_mixes[0]  # assume only 1 layer's dist
            for mix in route:
                p = probs[mix_list.index(mix)]
                if p > 0:
                    # entropy -= p * np.log2(p)
                    entropy += -np.log2(p + epsilon) + np.random.uniform(0, 0.01)

        return entropy

    def softmax(self, x):
        e_x = np.exp(np.array(x) - np.max(x))
        return e_x / e_x.sum()