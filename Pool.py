from random import sample
from Client import Client
from Mix import Mix

flushed = False  # True if a mix has flushed since the beginning of the simulation


class Pool(Mix):

    def __init__(self, mix_id, simulation, position,threshold, flush_percent,n_targets, corrupt, pr_mix):
        super().__init__(mix_id, simulation, position, n_targets, corrupt)

        self.pool = []
        self.pr_mix= pr_mix
        self.neighbors = []  # mixes in next layer used if routing == 'hopbyhop'
        self.threshold = threshold  # pool threshold
        self.flush_percent = flush_percent
        self.round = 0

    def receive_message(self, msg):
        """
        Handle both transit messages and destination messages
        """
        # Check if this mix is the final destination (only when mix_as_client is enabled)
        if hasattr(self.simulation, 'mix_as_client') and self.simulation.mix_as_client:
            if msg.route[-1] == self:
                # This mix is the destination - act like a client
                print(f"[GRID MIX-DEST] Mix {self.id} at ({self.grid_row},{self.grid_col}) received message {msg.id} as destination")
                msg.timeReceived = self.env.now
                self.log.received_messages_f(msg)
                
                if msg.target_bool and self.simulation.printing:
                    print(f'[GRID MIX-DEST] Target message {msg.id} arrived at destination mix {self.id}')
                
                # Send acknowledgment back to sender
                if hasattr(msg.route[0], 'receive_ack'):
                    msg.route[0].receive_ack(msg)
                return
        msg.next_hop_index += 1
        self.pool.append(msg)
        for i in range(0, self.n_targets):
            self.Pmix[i] += msg.pr_target[i]
        if msg.target_bool and self.simulation.printing:
            print(
                f'Target message arrived at mix {self.id} at time {self.env.now} and size of the pool{len(self.pool)}')
        if len(self.pool) >= self.threshold:
            print(f"[DEBUG] Threshold reached at Mix {self.id}")
            self.flush()
            self.round +=1
            self.simulation.numberrounds.append(self.round)
            for i in range(self.simulation.n_mixes_per_layer):
                self.simulation.set_stable_mix(i)

    def flush(self):
        flush_amount = int(self.flush_percent * self.threshold)
        flushing_list = sample(self.pool, k=flush_amount)

        for message in flushing_list:
            print(f"[DEBUG] Flushing Msg {message.id} at Mix {self.id}")
            self.messages_processed += 1
            self.update_probabilities(message)
            if not isinstance(message.route[message.next_hop_index], Client) and message.route[message.next_hop_index] is None:
                # not the last mix, for hop by hop routing
                message.route[message.next_hop_index] = sample(self.neighbors, k=1)[0]

            next_hop_index = message.route[message.next_hop_index]
            self.pool.remove(message)
            self.env.process(self.simulation.attacker.relay(message, next_hop_index))

    def update_probabilities(self, msg):
        if not self.corrupt:
            for j in range(0, self.n_targets):
                msg.pr_target[j] = self.Pmix[j] / len(self.pool)
                self.Pmix[j] = self.Pmix[j] - msg.pr_target[j]
    
    def setup_as_client(self, simulation, rate_client, mu, probability_dist_mixes, n_targets, n_hops, client_dummies, rate_client_dummies, Log):
        """
        Set up this pool mix to also act as a client (send messages) - only if enabled
        """
        # Client-related attributes
        self.rate_client = rate_client
        self.mu_client = mu  # Rename to avoid conflict with mix mu
        self.probability_dist_mixes = probability_dist_mixes
        self.n_targets = n_targets
        self.n_hops = n_hops
        self.client_dummies = client_dummies
        self.rate_client_dummies = rate_client_dummies
        self.log = Log
        self.message_id = 1
        self.other_clients = set()  # Will be set later
        
        # Get all other grid mixes for routing
        self.all_mixes = []
        if simulation.topology == 'grid':
            for layer in range(1, len(simulation.network.network_dict) + 1):
                self.all_mixes += simulation.network.network_dict[layer]
            # Remove self from available targets
            self.all_mixes = [mix for mix in self.all_mixes if mix.id != self.id]
        
        # Start message sending process
        self.env.process(self.send_message('Real', self.rate_client))
        print(f"[GRID MIX-CLIENT] Mix {self.id} at ({self.grid_row},{self.grid_col}) started sending messages")
        
        if self.client_dummies:
            self.env.process(self.send_message('ClientDummy', self.rate_client_dummies))
            print(f"[GRID MIX-CLIENT] Mix {self.id} started dummy message sending")

    def send_message(self, message_type, rate_client):
        """
        Send messages like a client (only active when mix_as_client is enabled)
        """
        from numpy.random import exponential
        
        while True:
            message, delay = self.create_message_as_sender(message_type, rate_client)
            yield self.env.timeout(delay)
            self.messages_sent += 1
            print(f"[GRID MIX-CLIENT] Mix {self.id} sending message {message.id} after delay {delay}")
            message.time_left = self.env.now
            self.log.sent_messages_f(message)
            
            # Send to next hop
            self.env.process(self.simulation.attacker.relay(message, message.route[1]))

    def create_message_as_sender(self, message_type, rate_client):
        """
        Create message with this mix as the sender (grid routing logic)
        """
        from numpy.random import exponential
        from random import choice
        import numpy as np
        
        np.random.seed()
        delay_client = exponential(rate_client)
        route = [self]  # Start with this mix as sender
        route_ids = [self.id]
        delays = [delay_client]
        link_delays = []
        pr_target = [0.0 for _ in range(self.n_targets)]
        
        # Grid topology routing (same logic as Client.py)
        current_node = self  # Start from this mix
        
        print(f"[GRID MIX-CLIENT] Mix {self.id} at ({self.grid_row},{self.grid_col}) creating route")
        
        # Create path through grid
        for hop in range(self.n_hops - 1):
            if not current_node.neighbors:
                break
                
            # Get unvisited neighbors
            unvisited_neighbors = [neighbor for neighbor in current_node.neighbors 
                                if neighbor not in route]
            
            if unvisited_neighbors:
                next_node = choice(unvisited_neighbors)
                print(f"[GRID MIX-CLIENT] Hop {hop+2}: {current_node.id}({current_node.grid_row},{current_node.grid_col}) -> {next_node.id}({next_node.grid_row},{next_node.grid_col}) [unvisited]")
            else:
                # Allow backtracking but avoid immediate previous
                if len(route) >= 2:
                    available_neighbors = [neighbor for neighbor in current_node.neighbors 
                                        if neighbor != route[-2]]
                else:
                    available_neighbors = current_node.neighbors
                
                if available_neighbors:
                    next_node = choice(available_neighbors)
                    print(f"[GRID MIX-CLIENT] Hop {hop+2}: {current_node.id} -> {next_node.id} [backtracking]")
                else:
                    next_node = choice(current_node.neighbors)
                    print(f"[GRID MIX-CLIENT] Hop {hop+2}: {current_node.id} -> {next_node.id} [forced]")
            
            delay_per_mix = exponential(self.mu_client)
            delays.append(delay_per_mix)
            route.append(next_node)
            route_ids.append(next_node.id)
            current_node = next_node
        
        # Add final delay and choose destination mix
        delays += [0]
        destination_mix = choice(list(self.other_clients))  # Choose another mix as destination
        route += [destination_mix]
        route_ids += [destination_mix.id]
        
        # Build link delays
        for i in range(len(route) - 1):
            link_delay = 0.1
            link_delays.append(link_delay)
            print(f"[GRID MIX-CLIENT] Link Delay Hop {i}: {route[i].id} -> {route[i+1].id}: {link_delay}")
        
        from Message import Message
        message = Message(self.message_id, message_type, self, route, delays, link_delays, pr_target, False)
        
        print(f"[GRID MIX-CLIENT] Mix {self.id} created message {message.id}")
        print(f"[GRID MIX-CLIENT] Route: {[node.id for node in route]}")
        
        self.message_id += 1
        return message, delay_client

    def receive_ack(self, message):
        """Handle acknowledgment (from Client interface)"""
        print(f"[GRID MIX-CLIENT] Mix {self.id} received ACK for message {message.id}")

    def client_str(self):
        """String representation when acting as client"""
        return f'Client id: {self.id}'