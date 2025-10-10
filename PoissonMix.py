from random import sample, choice
from Client import Client
from Mix import Mix
from numpy.random import exponential
from Message import Message
import numpy as np
import random

global idDummies


class PoissonMix(Mix):
    def __init__(self, mix_id, simulation, position,link_based_dummies,multiple_hop_dummies, rate_mix_dummies, n_targets, corrupt, pr_mix):
        super().__init__(mix_id, simulation, position, n_targets, corrupt)
        self.pool = []  # Pool
        self.neighbors = set()
        self.link_based_dummies = link_based_dummies
        self.multiple_hop_dummies = multiple_hop_dummies
        self.rate_mix_dummies = rate_mix_dummies
        self.pr_mix = pr_mix
        self.pool_dummies = []
        if (self.link_based_dummies or self.multiple_hop_dummies) and (not self.corrupt) and self.layer != self.simulation.n_layers:
            self.env.process(self.send_dummies())
        # Add layer position mapping
        # self.layer_position = self.get_layer_position()
    
    #for larmix stratified topology
    def get_layer_position(self):
        """Get the position index of this mix within its layer (0-based)"""
        if self.simulation.topology == 'stratified':
            # Find this mix's position in its layer
            layer_mixes = self.simulation.network.network_dict[self.layer]
            for i, mix in enumerate(layer_mixes):
                if mix.id == self.id:
                    return i
        return 0  # fallback

    def receive_message(self, msg):
        if not self.simulation.startAttack:  # if a mix reaches a poolsize of 5 percent higher than the average,
            # the GPA can monitor the network and choose a target message
            clients = self.simulation.n_clients
            lambdaClient = self.simulation.rate_client
            average = (clients * lambdaClient * self.simulation.mu) * self.pr_mix
            var1 = len(self.pool) >= average
            print(f"var1 -> {var1}  len(self.pool) -> {len(self.pool)}  average -> {average}")
            if var1:
                if self.simulation.topology == 'stratified':
                    # print(f'[{msg.id}] Entered mix {self.id} ReceiveMsg() at time {self.env.now}')
                    if self.layer == 1:
                        if self.simulation.routing == 'larmix':
                            layer_position = self.get_layer_position()
                            print(f'[{msg.id}] [Mix {self.id}] Calling set_stable_mix for index {layer_position} at time {self.env.now}')
                            self.env.process(self.simulation.set_stable_mix(layer_position))
                        else:
                            self.env.process(self.simulation.set_stable_mix(self.id - 1))
                        
                    #if all(self.simulation.stableMixL1):
                        #for i in range(len(self.simulation.stableMixL1)):
                            #self.simulation.setStableMix(i)
                elif self.simulation.topology == 'cyclic_stratified':
                    if not self.simulation.stable_layer[self.layer - 1]:
                        # mark this layer stable *once*
                        self.simulation.stable_layer[self.layer - 1] = True
                        if self.simulation.printing:
                            print(f"[{self.env.now}] Layer {self.layer} stable "
                                f"({sum(self.simulation.stable_layer)}/"
                                f"{self.simulation.n_layers})")
                        if all(self.simulation.stable_layer):
                            self.simulation.startAttack = True
                            if self.simulation.printing:
                                print(f"[{self.env.now}] Ring stable → startAttack = True")
                elif self.simulation.topology == 'XRD':
                    if var1:
                        self.env.process(self.simulation.setStableChain(self.n_chain))
                    if all(self.simulation.stableChains):
                        for i in range(len(self.simulation.stableChains)):
                            self.simulation.setStableChain(i)
                elif self.simulation.topology == 'free route':
                    if self.simulation.routing == 'larmix':
                        # For LARMix free route, use mix position in the all_mixes list
                        all_mixes = self.simulation.network.network_dict[1]
                        mix_position = next((i for i, mix in enumerate(all_mixes) if mix.id == self.id), 0)
                        print(f"[{self.env.now}] LARMix Free Route: Mix {self.id} (position {mix_position}) stable")
                        self.env.process(self.simulation.set_stable_mix(mix_position))
                    else:
                        # Original logic for non-LARMix free route
                        print(f"[{self.env.now}] Mix {self.id} stable => set_stable_mix({self.id - 1})")
                        self.env.process(self.simulation.set_stable_mix(self.id - 1))
            if self.simulation.topology == 'ba topology':
                if len(self.pool) >= 2:  
                    print(f"[BA Debug] Mix {self.id} is stable => Calling set_stable_mix({self.id - 1})")
                    self.env.process(self.simulation.set_stable_mix(self.id - 1))
        for i in range(0, self.n_targets):
            self.Pmix[i] += msg.pr_target[i]
        if msg.target_bool and self.simulation.printing:
            print(f'[{msg.id}] [Mix {self.id}] Target message arrived at time {self.env.now} and Number of messages inside '
                  f'the pool {len(self.pool)}')
        msg.next_hop_index += 1
        if msg.route[msg.next_hop_index] == None:
            msg.route[msg.next_hop_index] = random.choice(list(self.neighbors))
        if msg.type == 'Real':
            self.pool.append(msg)
            self.env.process(self.send_msg(msg))
        elif msg.type == 'Dummy':
            if self.link_based_dummies:
                self.drop_dummies(msg)
            elif self.multiple_hop_dummies:
                if self.layer == self.simulation.n_layers:
                    self.drop_dummies(msg)
                elif self.layer != self.simulation.n_layers:
                    self.pool.append(msg)
                    self.env.process(self.send_msg(msg))

    def send_msg(self, msg):
        hop_index = msg.next_hop_index - 1  # because next_hop_index starts at 1 after leaving client
        if hop_index < 0:
            hop_index = 0  # safety for first hop
        delay_for_this_hop = msg.delays[hop_index] 
        
        yield self.env.timeout(delay_for_this_hop)
        print(f"[{msg.id}] [Hop {hop_index}] [Mix {self.id}] Processing Delay Completed => {delay_for_this_hop}")
        self.messages_processed += 1
        self.update_probabilities(msg, len(self.pool))
        next_hop = msg.route[msg.next_hop_index]
        self.pool.remove(msg)
        self.env.process(self.simulation.attacker.relay(msg, next_hop))

    def update_probabilities(self, msg, pool_size):
        if not self.corrupt:
            for j in range(0, self.n_targets):
                msg.pr_target[j] = self.Pmix[j] / pool_size
                self.Pmix[j] = self.Pmix[j] - msg.pr_target[j]

    def drop_dummies(self, msg):
        if self.layer == self.simulation.n_layers:
            self.simulation.Log.dummies_dropped_end_link( msg, self.id)

    def send_dummies(self):
        dummy_id = 1
        while True:
            new_message = self.create_dummies(dummy_id)
            new_message.creator = self.id
            self.pool.append(new_message)
            yield self.env.timeout(exponential(self.rate_mix_dummies))
            dummy_id += 1
            self.env.process(self.send_msg(new_message))