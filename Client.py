from random import choice, sample
from Message import Message
from numpy.random import exponential
import numpy as np
import random
import itertools
from BatchTracker import (incoming_batches, 
                          next_incoming_batch_id, 
                          outgoing_batches, 
                          next_outgoing_batch_id, 
                          incoming_outgoing_batch_map, 
                          outgoing_to_incoming_batch_map,
                          compute_batch_permutations,
                          )

class Client:
    def __init__(self, simulation, id, network_dict, rate_client, mu, probability_dist_mixes, n_targets, n_hops, client_dummies, rate_client_dummies, Log, batch_size):
        self.id = id
        self.env = simulation.env
        self.simulation = simulation  # simulation object
        self.network_dict = network_dict
        self.class_ends = self.env.event()
        self.mu = mu  # Avg delay: for delays at the poisson mixes
        self.probability_dist_mixes = probability_dist_mixes
        self.other_clients = set()
        self.message_id = 1
        #batch-algorithm
        self.current_batch_id = None
        self.sent_msg_count_in_batch = 0
        self.current_batch_receiver = None

        self.rate_client = rate_client
        self.all_mixes = []
        self.n_targets = n_targets
        self.n_hops = n_hops
        self.batch_size = batch_size
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
        assert self.current_batch_receiver is not None, "Receiver for batch not set!"
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

        elif self.simulation.topology == 'cyclic_stratified':
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
        # batch algorithm
        receiver = self.current_batch_receiver
        # receiver = sample(list(self.other_clients), k=1)[0]
        print(f"[Route Delays]: {delays}")
        print(f"==>> Receiver: {receiver} at time {self.env.now}")
        route += [receiver]
        route_ids += [receiver.id]
        print(f"[Debug] ==>> Route: {route}") 

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
        global next_outgoing_batch_id
        message.timeReceived = self.env.now

        # batch algorithm
        incoming_batch_id = message.incoming_batch_id

        if incoming_batch_id not in incoming_outgoing_batch_map:
            out_batch_id = next_outgoing_batch_id
            next_outgoing_batch_id += 1
            incoming_outgoing_batch_map[incoming_batch_id] = out_batch_id
            print(f"==>> Mapping IncBatch {incoming_batch_id} to OutBatch {out_batch_id}")
            if out_batch_id not in outgoing_to_incoming_batch_map:
                outgoing_to_incoming_batch_map[out_batch_id] = incoming_batch_id
                print(f"==>> Mapping OutBatch {out_batch_id} to IncBatch {incoming_batch_id}")
        else:
            out_batch_id = incoming_outgoing_batch_map[incoming_batch_id]
        # print(f"==>> Inc to Out Batch Map: {incoming_outgoing_batch_map}")

        # Extract incoming msg number from msg id (format: M_batchid_msgno)
        incoming_msg_id = message.incoming_msg_id
        incoming_msg_no = incoming_msg_id.split('_')[-1]

        # Assign outgoing msg id as O_outbatchid_incomingmsgno
        out_msg_id = f"O_{out_batch_id}_{incoming_msg_no}"
        message.outgoing_batch_id = out_batch_id
        message.outgoing_msg_id = out_msg_id
        print(f'IncomingMsgID: {incoming_msg_id}\nOutgoingMsgID: {out_msg_id}')

        # Update global outgoing_batches dict
        if out_batch_id not in outgoing_batches:
            outgoing_batches[out_batch_id] = {}
        outgoing_batches[out_batch_id][out_msg_id] = message.timeReceived
        print(f"==>> {out_msg_id} Received at : {message.timeReceived}")
        # print(f"==>> Outgoing Batches: {outgoing_batches}")
        
        self.log.received_messages_f(message)
        if message.target_bool and self.simulation.printing:
            print(f'Target message arrived at destination Client at time {self.env.now}')
        if message.type == 'Real' or message.type == 'ClientDummy':
            message.route[0].receive_ack(message)
        # Compute and print all possible permutations for each outgoing batch
        compute_batch_permutations(self, message)

    def send_message(self, message_type, rate_client):
        global next_incoming_batch_id
        while True:
            # batch-algorithm
            # If not currently sending a batch, claim the next available batch id
            if self.current_batch_id is None or self.sent_msg_count_in_batch >= self.batch_size:
                self.current_batch_id = next_incoming_batch_id
                next_incoming_batch_id += 1
                self.sent_msg_count_in_batch = 0
                self.current_batch_receiver = sample(list(self.other_clients), k=1)[0]
                
            batch_id = self.current_batch_id
            msg_number = self.sent_msg_count_in_batch
            msg_id = f"M_{batch_id}_{msg_number}"
            

            message, delay = self.create_message(message_type, rate_client)
            message.incoming_batch_id = batch_id
            message.incoming_msg_id = msg_id

            yield self.env.timeout(delay)
            print(f"==>> Incoming Msg id: {msg_id}")
            print(f"==>> Sending Delay: {delay}")
            message.time_left = self.env.now

            # Track in global incoming_batches
            if batch_id not in incoming_batches:
                incoming_batches[batch_id] = {}
            incoming_batches[batch_id][msg_id] = message.time_left
            print(f"==>> {msg_id} Left at : {message.time_left}")
            print(f"==>> Incoming Batches: {incoming_batches}")

            self.log.sent_messages_f(message)
            self.env.process(self.simulation.attacker.relay(message, message.route[1]))
            print(f"==>> message.route[1]: {message.route}")

            # Update per-client batch message count
            self.sent_msg_count_in_batch += 1

            # If finished this batch, next loop will claim a new batch id
            if self.sent_msg_count_in_batch >= self.batch_size:
                self.current_batch_id = None

    def receive_ack(self, message):  # Message received
        pass

    def __str__(self):
        return 'Client id: {}'.format(self.id)

    def __repr__(self):
        return self.__str__()

    

        