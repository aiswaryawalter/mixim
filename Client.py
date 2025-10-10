from random import choice, sample
from Message import Message
from numpy.random import exponential
import numpy as np
import random

class Client:
    def __init__(self, simulation, id, network_dict, rate_client, mu, probability_dist_mixes, n_targets, n_hops, client_dummies, rate_client_dummies, Log):
        self.id = id
        self.server_id = None 
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
        # Link delays (between nodes) - will be filled based on route
        link_delays = []
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
            all_server_ids = [int(sid) for sid in self.simulation.node_coords.keys()]

            # Get server IDs already used by mix nodes
            mix_server_ids = set()
            for layer in range(1, self.simulation.network.num_layers + 1):
                for mix in self.simulation.network.network_dict[layer]:
                    mix_server_ids.add(int(mix.server_id))
            
            # print(f"[DEBUG] All server IDs: {sorted(all_server_ids)}")
            # print(f"[DEBUG] Mix server IDs: {sorted(mix_server_ids)}")
            
            # Available server IDs for clients (excluding mix servers)
            available_client_server_ids = [int(sid) for sid in all_server_ids if sid not in mix_server_ids]
            
            # print(f"[DEBUG] Available client server IDs: {sorted(available_client_server_ids)}")
            
            if len(available_client_server_ids) < 2:
                raise ValueError(f"Not enough available server IDs for clients. Need at least 2, have {len(available_client_server_ids)}")

            sender_server_id = np.random.choice(available_client_server_ids) 
            self.server_id = sender_server_id 
            receiver_server_candidates = [sid for sid in available_client_server_ids if sid != sender_server_id]
            receiver_server_id = np.random.choice(receiver_server_candidates)
            receiver = sample(list(self.other_clients), k=1)[0]
            receiver.server_id = receiver_server_id

            print(f"[LARMix] Sender Client: {self.id} --- Server: {sender_server_id}")
            print(f"[LARMix] Receiver Client: {receiver.id} --- Server: {receiver_server_id}")

            best_route = self.simulation.network.larmix.sample_larmix_path_stratified(self.simulation.tau, sender_server_id, receiver_server_id)
            route = [self] + best_route + [receiver]
            route_ids = [node.id for node in route]
             
            # Processing delays at each mix
            for _ in range(len(best_route)):
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
            
            print(f"[LARMix] Complete Route: {[(getattr(node, 'server_id', node.server_id), getattr(node, 'id', node.id)) for node in route]}")

        elif (self.simulation.routing == 'larmix'
              and self.simulation.topology == 'free route'):
            # NEW: LARMix with free route topology
            all_server_ids = [int(sid) for sid in self.simulation.node_coords.keys()]

            # Get server IDs used by mix nodes (all mixes are in network_dict[1] for free route)
            mix_server_ids = set()
            for mix in self.all_mixes:  # all_mixes contains all mixes in free route
                if hasattr(mix, 'server_id'):
                    mix_server_ids.add(int(mix.server_id))
            
            # Available server IDs for clients
            available_client_server_ids = [int(sid) for sid in all_server_ids if sid not in mix_server_ids]
            
            if len(available_client_server_ids) < 2:
                raise ValueError(f"Not enough available server IDs for clients. Need at least 2, have {len(available_client_server_ids)}")

            # Assign server IDs to sender and receiver
            sender_server_id = np.random.choice(available_client_server_ids) 
            self.server_id = sender_server_id 
            receiver_server_candidates = [sid for sid in available_client_server_ids if sid != sender_server_id]
            receiver_server_id = np.random.choice(receiver_server_candidates)
            receiver = sample(list(self.other_clients), k=1)[0]
            receiver.server_id = receiver_server_id

            print(f"[LARMix Free Route] Sender Client: {self.id} --- Server: {sender_server_id}")
            print(f"[LARMix Free Route] Receiver Client: {receiver.id} --- Server: {receiver_server_id}")

            # Sample latency-aware path through free route network
            best_route = self.simulation.network.larmix.sample_larmix_path_free_route(
                self.simulation.tau, sender_server_id, receiver_server_id, self.n_hops
            )
            route = [self] + best_route + [receiver]
            route_ids = [getattr(node, 'id', node.id) for node in route]
            
            # Processing delays at each mix
            for _ in range(len(best_route)):
                delay_per_mix = exponential(self.mu)
                delays.append(delay_per_mix)
            
            print(f"[LARMix Free Route] Route: {[(getattr(node, 'server_id', 'Client'), getattr(node, 'id', node.id)) for node in route]}")

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
        if self.simulation.routing != 'larmix':
            receiver = sample(list(self.other_clients), k=1)[0]
            # print(f"[Processing Delays]: {delays}")
            print(f"==>> Receiver chosen: {receiver} at time {self.env.now}")
            route += [receiver]
            route_ids += [receiver.id]
            print(f"[Route]: {route}") 

        # Build link delays based on route
        if self.simulation.routing == 'larmix':
            for i in range(len(route) - 1):
                prev_node = route[i]
                next_node = route[i + 1]
            
                prev_server_id = getattr(prev_node, 'server_id', None)
                next_server_id = getattr(next_node, 'server_id', None)
                
                if prev_server_id and next_server_id:
                    # Use latency matrix for LARMix
                    link_delay = self.simulation.network.latency_matrix.get(
                        (prev_server_id, next_server_id), 0.05
                    )
                    # print(f"[Link Delay] ServerID {prev_server_id} -> {next_server_id}: {link_delay}")
                else:
                    # Default link delay
                    link_delay = 0.05
                
                link_delays.append(link_delay)
                # print(f"[Link Delay] ServerID {prev_server_id} -> {next_server_id}: {link_delay}")
        else:
            # All other routing strategies use fixed link delay of 0.05
            for i in range(len(route) - 1):
                link_delay = 0.1 #avg of latencies.csv
                link_delays.append(link_delay)
                print(f"[Link Delay] Hop {i}: {route[i].id} -> {route[i+1].id}: {link_delay}")
        
            
        message = Message(self.message_id, message_type, self, route, delays, link_delays, pr_target,False)
        if self.message_id == 1 and self.id ==1:
            for i in range(len(self.probability_dist_mixes)):
                if self.simulation.printing:
                    print("Weights Layer %d %s"%(i, self.probability_dist_mixes))
                else:
                    pass
        print(f"[{message.id}] [Link Delays]: {link_delays}")
        print(f"[{message.id}] [Processing Delays]: {delays}")
        print(f"[{message.id}] [Route]: {route}")
        self.message_id += 1
        return message, delay_client



    def receive_message(self, message):
        message.timeReceived = self.env.now
        self.log.received_messages_f(message)
        if message.target_bool and self.simulation.printing:
            print(f'[{message.id}] [Hop {message.next_hop_index}] [Client {self.id}] Target message arrived at destination at time {self.env.now}')
        if message.type == 'Real' or message.type == 'ClientDummy':
            message.route[0].receive_ack(message)


    def send_message(self, message_type, rate_client):
        while True:
            message, delay = self.create_message(message_type, rate_client)
            yield self.env.timeout(delay)
            print(f"==>> [{message.id}] [Client {self.id}] Processing Delay at Sender Client Completed : {delay}")
            message.time_left = self.env.now
            self.log.sent_messages_f(message)

            self.env.process(self.simulation.attacker.relay(message, message.route[1]))
            print(f"==>> [{message.id}] Relaying to the next mix node: {message.route[2]}")

    def receive_ack(self, message):  # Message received
        pass

    def __str__(self):
        return 'Client id: {}'.format(self.id)

    def __repr__(self):
        return self.__str__()
    