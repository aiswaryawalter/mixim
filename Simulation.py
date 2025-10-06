import simpy
import Client
from Network import Network
import pandas as pd
import numpy as np
from Relay import Attacker
from Log import Log
from util import XRD_New, load_latency_matrix, load_server_coords
import os

DEFAULT_TOPOLOGY = 'stratified'
logDir = 'Logs/'
if not os.path.exists(logDir):
    os.makedirs(logDir)


class Simulation(object):

    def __init__(self, mix_type, simDuration, rate_client, mu, logging, topology, fully_connected, n_clients, n_hops, 
                 flush_percent, printing, flush_timeout, threshold, routing, latency_bound, tau, balancing,
                 n_layers, n_mixes_per_layer, corrupt, unifrom_corruption, probability_dist_mixes, nbr_cascacdes, m_barabasi_mixes, client_dummies,
                 rate_client_dummies, link_based_dummies, multiple_hops_dummies, rate_mix_dummies, Network_template):

        self.Log = Log()
        self.logs = []
        self.logging = logging
        self.printing = printing
        self.topology = topology
        self.n_cascades = nbr_cascacdes
        self.fully_connected = fully_connected
        self.flush_percent = flush_percent

        self.client_dummies = client_dummies
        self.rate_client_dummies = rate_client_dummies
        self.link_based_dummies = link_based_dummies
        self.multiple_hop_dummies = multiple_hops_dummies
        self.rate_mix_dummies = rate_mix_dummies

        #  larmix
        self.latency_matrix = load_latency_matrix("latency.csv")
        self.node_coords = load_server_coords("servers.csv")
        self.tau = tau
        self.balancing = balancing
        self.latency_bound = latency_bound

        self.n_clients = n_clients
        self.n_hops = n_hops
        self.clientsSet = set()
        self.rate_client = rate_client  # average delay between messages being sent from client
        self.threshold = threshold
        self.mu = mu  # average delay at poisson mixes
        self.n_layers = n_layers
        self.n_mixes_per_layer = n_mixes_per_layer
        self.m_barabasi_mixes = m_barabasi_mixes
        self.corrupt = corrupt
        self.probability_dist_mixes = probability_dist_mixes
        self.unifrom_corruption = unifrom_corruption
        self.mix_type = mix_type
        self.routing = routing
        self.env = simpy.Environment()
        self.SimDuration = simDuration
        self.burnout = 10
        self.flush_timeout = flush_timeout
        self.n_targets = 0
        self.MsgsDropped = []

        

        self.dummyID = 0
        time_stable = ((1 / self.rate_client) / self.n_layers) * self.mu + 2
        if self.mix_type == 'poisson':
            self.n_targets = int(((self.SimDuration - time_stable) ) / 2)
            print(f"n_targets={self.n_targets}")
        else:
            self.n_targets = int((self.SimDuration - self.flush_timeout - 1) / 4)
        self.network = Network(self.mix_type, self.n_layers, self.n_mixes_per_layer, 
                               self.corrupt,self.unifrom_corruption, self, self.threshold,
                               self.flush_percent, self.topology, fully_connected, self.flush_timeout,
                               self.probability_dist_mixes,
                               self.n_cascades, self.m_barabasi_mixes, self.link_based_dummies, self.multiple_hop_dummies,
                               self.rate_mix_dummies,
                               Network_template, self.n_targets)

        self.set_clients(self.probability_dist_mixes, self.n_targets, self.client_dummies, self.rate_client_dummies,
                         self.Log)
        # self.stableMix = [False for i in range(self.n_mixes_per_layer*self.n_layers)]  # only start attack after mixes are stable
        self.stableChains = [False for i in range(1, 1 + 6)]  # only start attack after chains are stable
        if self.topology == 'stratified':
            self.stableMixL1 = [False for i in range(self.n_mixes_per_layer)]  # only start attack after mixes are stable
        if self.topology == 'cyclic_stratified':
            self.stable_layer = [False] * self.n_layers 
        self.attacker = Attacker(self, self.n_targets)  # attacker/relay object
        self.endEvent = self.env.event()  # event that triggers the end of the simulation
        self.TargetMessageEnd = False  # if target message has reached the end client
        self.startAttack = False  # if the attacker is allowed to choose a target message
        self.NumberMsgsDropped = 0
        self.numberrounds = []
        # Start periodic log saving process
        # self.env.process(self.periodic_log_saver())

    def set_stable_mix(self, index):
        print(f"[{self.env.now}] Entered set_stable_mix for index={index}")
        if self.mix_type == 'pool':
            yield self.env.timeout(10)
            self.startAttack = True
        elif self.mix_type == 'time':
            yield self.env.timeout(self.flush_timeout + 5)
            self.startAttack = True
        self.stableMixL1[index] = True
        print(f"[{self.env.now}] stableMixL1[{index}] = {self.stableMixL1[index]}")
        if all(self.stableMixL1):
            yield self.env.timeout(2)
            self.startAttack = True
        # if self.routing == 'larmix':
        #     yield self.env.timeout(2)
        #     self.startAttack = True
        print(f"[{self.env.now}] set_stable_mix done => startAttack={self.startAttack}")


    def set_stable_chain(self, position):
        if self.mix_type == 'pool':
            yield self.env.timeout(10)
            self.startAttack = True
        self.stableChains[position - 1] = True
        if all(self.stableChains):
            yield self.env.timeout(2)
            self.startAttack = True

    def set_clients(self, probabilityDistribution, n_targets, client_dummies, rate_client_dummies, Log):
        if self.topology == 'stratified':
            for client_no in range(self.n_clients):
                client = Client.Client(self, client_no, self.network.network_dict, self.rate_client, self.mu,
                                       probabilityDistribution, n_targets, self.n_hops, client_dummies, rate_client_dummies, Log)
                self.clientsSet.add(client)
            for client in self.clientsSet:
                client.other_clients = self.clientsSet - {client}
        
        elif self.topology == 'cyclic_stratified':
            for client_no in range(self.n_clients):
                client = Client.Client(
                    self,
                    client_no,
                    self.network.network_dict,   # the ring
                    self.rate_client,
                    self.mu,
                    probabilityDistribution,
                    n_targets,
                    self.n_hops,
                    client_dummies,
                    rate_client_dummies,
                    Log
                )
                self.clientsSet.add(client)

            for client in self.clientsSet:
                client.other_clients = self.clientsSet - {client}

        elif self.topology == 'XRD':
            groups_lists = XRD_New(self.network.list_cascades)
            n_group_client = self.n_clients // len(groups_lists)
            for client_id in range(n_group_client):
                client = Client.Client(self, client_id, groups_lists[0], self.rate_client, self.mu,
                                       probabilityDistribution, n_targets, self.n_hops, client_dummies, Log)
                self.clientsSet.add(client)
            for n_client in range(n_group_client, n_group_client * 2):
                client = Client.Client(self, n_client, groups_lists[1], self.rate_client, self.mu,
                                       probabilityDistribution, n_targets, self.n_hops, client_dummies, Log)
                self.clientsSet.add(client)
            for n_client in range(n_group_client * 2, n_group_client * 3):
                client = Client.Client(self, n_client, groups_lists[2], self.rate_client, self.mu,
                                       probabilityDistribution, n_targets, self.n_hops, client_dummies, Log)
                self.clientsSet.add(client)
            for n_client in range(n_group_client * 3, self.n_clients):
                client = Client.Client(self, n_client, groups_lists[3], self.rate_client, self.mu,
                                       probabilityDistribution, n_targets, self.n_hops, client_dummies, Log)
                self.clientsSet.add(client)
            for client in self.clientsSet:
                client.otherClients = self.clientsSet - {client}

        elif self.topology == 'free route':
            for client_no in range(self.n_clients):
                client = Client.Client(
                    self,
                    client_no,
                    self.network.network_dict,  # So they can see all mixes
                    self.rate_client,
                    self.mu,
                    probabilityDistribution,
                    n_targets,
                    self.n_hops,
                    client_dummies,
                    rate_client_dummies,
                    Log
                )
                self.clientsSet.add(client)
            for client in self.clientsSet:
                client.other_clients = self.clientsSet - {client}

        elif self.topology == 'ba topology':
            for client_no in range(self.n_clients):
                client = Client.Client(
                    self,
                    client_no,
                    self.network.network_dict, 
                    self.rate_client,
                    self.mu,
                    probabilityDistribution,
                    n_targets,
                    self.n_hops,
                    client_dummies,
                    rate_client_dummies,
                    Log
                )
                self.clientsSet.add(client)
            for client in self.clientsSet:
                client.other_clients = self.clientsSet - {client}

    def periodic_log_saver(self):
        """Periodically save logs during simulation"""
        while True:
            try:
                yield self.env.timeout(100)  # Save every 100 time units
                current_time = self.env.now
                
                # Only save if there's new data
                if (len(self.Log.sent_messages["MessageID"]) > 0 or 
                    len(self.Log.received_messages["MessageID"]) > 0):
                    
                    print(f"[{current_time}] Periodic log save...")
                    self.save_current_logs()
                    
            except Exception as e:
                print(f"[ERROR] Periodic logging failed: {e}")
                break

    def save_current_logs(self):
        """Save current state of logs"""
        try:
            df_sent_messages = pd.DataFrame(self.Log.sent_messages)
            df_received_messages = pd.DataFrame(self.Log.received_messages)
            df_dummies_messages = pd.DataFrame(self.Log.dummy_messages)
            
            timestamp = int(self.env.now)
            
            print(f"[LOG DEBUG] Current simulation time: {timestamp}")
            print(f"[LOG DEBUG] Sent messages count: {len(df_sent_messages)}")
            print(f"[LOG DEBUG] Received messages count: {len(df_received_messages)}")
            print(f"[LOG DEBUG] Dummy messages count: {len(df_dummies_messages)}")
            
            # Save even if empty (with headers)
            df_sent_messages.to_csv(f'{logDir}SentMessages_t{timestamp}.csv', index=False)
            print(f"[LOG DEBUG] Saved sent messages to SentMessages_t{timestamp}.csv")
            
            df_received_messages.to_csv(f'{logDir}ReceivedMessages_t{timestamp}.csv', index=False)
            print(f"[LOG DEBUG] Saved received messages to ReceivedMessages_t{timestamp}.csv")
            
            df_dummies_messages.to_csv(f'{logDir}DummyMessages_t{timestamp}.csv', index=False)
            print(f"[LOG DEBUG] Saved dummy messages to DummyMessages_t{timestamp}.csv")
            
        except Exception as e:
            print(f"[LOG ERROR] Failed to save logs: {e}")
            import traceback
            traceback.print_exc()

    def run(self, time=None):
        if self.printing:
            print('\n')
            print('----------Simulation Data----------')
            print('Topology: {}'.format(self.topology))
            print('Routing strategy: {}'.format(self.routing))
            print('Mix type: {}'.format(self.mix_type))
            print('Layers: {}, amount of mixes per layer: {}'.format(self.n_layers, self.n_mixes_per_layer))
            print(
                'Amount of clients: {}, average delay between 2 messages: {}'.format(self.n_clients, self.rate_client))
            l = ''
            index = 0
            if self.topology == 'stratified':
                for layer in self.network.network_dict:
                    k = 'Layer {}: [ '.format(layer)
                    for mixnb in range(self.n_mixes_per_layer):
                        k += str(self.network.network_dict[layer][mixnb])
                        if self.network.network_dict[layer][mixnb].corrupt:
                            index += 1
                    l += k + ']'
                    l += '\n'
                print(l)
            elif self.topology == 'Cascade':
                k2 = 'Cascade'
                for cascade in self.network.list_cascades:
                    for mixnb in range(3):
                        k2 += str(cascade[mixnb].id)
                    l += k2 + ']'
                    l += '\n'
            print('----------Starting Simulation----------')
        if self.printing:
            print('Topology: {}'.format(self.topology))
        try:
            if time is None:
                self.env.run(until=self.endEvent)
            else:
                self.env.run(until=time)
        except KeyboardInterrupt:
            print(f"\n[INTERRUPTED][{self.env.now}] Simulation stopped by user. Saving current logs...")
            self.save_current_logs()
            print("Logs saved before exit.")
            return
        except Exception as e:
            print(f"\n[ERROR] Simulation failed: {e}. Saving current logs...")
            self.save_current_logs()
            raise
        finally:
            # Always save logs at the end, regardless of how simulation ended
            print("Saving final logs...")
            self.save_current_logs()

        if self.printing:
            print('----------Simulation Ended---------')
            print('\n')

        # Data from Clients(senders and receivers)
        df_sent_messages = pd.DataFrame(self.Log.sent_messages)
        df_received_messages = pd.DataFrame(self.Log.received_messages)
        df_dummies_messages = pd.DataFrame(self.Log.dummy_messages)

        if self.logging:
            df_sent_messages.to_csv(f'{logDir}SentMessages.csv')
            df_received_messages.to_csv(f'{logDir}ReceivedMessages.csv')
            df_dummies_messages.to_csv(f'{logDir}DummyMessages.csv')
        else:
            pass

        entropy = []
        for i in range(0, self.n_targets):
            entropy.append(0.0)
        tableProb = df_received_messages['MessageTarget'].to_numpy(copy=True)
        for j in range(0, self.n_targets):
            for m in range(len(tableProb)):
                if tableProb[m][j] != 0:
                    entropy[j] += - tableProb[m][j] * np.log2(tableProb[m][j])

        dict_entropy = {'Entropy': entropy}
        df_entropy = pd.DataFrame(dict_entropy)
        df_entropy.to_csv(f'{logDir}{self.n_layers}layers_{self.n_mixes_per_layer}mixes_player_Entropy.csv')

        entropy_mean = np.mean(entropy)
        try:
            entropy_median = np.median(entropy)
            entropy_q25 = np.quantile(entropy, .25)
        except:
            entropy_median = 0
            entropy_q25 = 0

        sum_delays = 0
        for re, le in zip(self.Log.received_messages["MessageTimeReceived"],
                          self.Log.received_messages["MessageTimeLeft"]):
            sum_delays += (re - le)
        average_delay = sum_delays / len(self.Log.received_messages["MessageTimeReceived"])

        # New latency calculations using the latency methods
        total_processing_latencies = []
        total_link_latencies = []
        total_latencies = []
        
        # Calculate latencies for each received message
        for processing_lat, link_lat, total_lat in zip(
            self.Log.received_messages["MessageProcessingLatency"],
            self.Log.received_messages["MessageLinkLatency"], 
            self.Log.received_messages["MessageTotalLatency"]):
            
            total_processing_latencies.append(processing_lat)
            total_link_latencies.append(link_lat)
            total_latencies.append(total_lat)
        
        # Calculate averages
        avg_processing_latency = sum(total_processing_latencies) / len(total_processing_latencies)
        avg_link_latency = sum(total_link_latencies) / len(total_link_latencies)
        avg_total_latency = sum(total_latencies) / len(total_latencies)
        
        # Calculate additional statistics
        max_processing_latency = max(total_processing_latencies)
        min_processing_latency = min(total_processing_latencies)
        max_link_latency = max(total_link_latencies)
        min_link_latency = min(total_link_latencies)
        max_total_latency = max(total_latencies)
        min_total_latency = min(total_latencies)
        
        if self.printing:
            print('----------Simulation Stats----------')
            print('\n')
            print('----------Simulation Data----------')
            print('Topology: {}'.format(self.topology))
            print('Routing strategy: {}'.format(self.routing))
            print('Mix type: {}'.format(self.mix_type))
            print('Mu: {}'.format(self.mu))
            print('Layers: {}, \n amount of mixes per layer: {}, \n Number of hops: {}'.format(self.n_layers, self.n_mixes_per_layer, self.n_hops))
            print(
                'Amount of clients: {}, \n average delay between 2 messages: {}'.format(self.n_clients, self.rate_client))
            # print('Average Latency: {}'.format(latency))
            print("Number of targets chosen", self.n_targets)
            print('Number of Real messages generated', len(self.Log.sent_messages["MessageID"]))
            print('Number of Real messages Received', len(self.Log.received_messages["MessageID"]))
            print('Number of Dummy messages dropped', len(self.Log.dummy_messages["DummyID"]))
            # print("Average delay per message", average_delay)
            # print('-------------------------------------')

            # Enhanced latency statistics
            print('\n----------Latency Statistics----------')
            print(f"Average time-based delay per message: {average_delay:.6f}")
            print(f"Average processing latency: {avg_processing_latency:.6f}")
            print(f"Average link latency: {avg_link_latency:.6f}")
            print(f"Average total latency: {avg_total_latency:.6f}")
            
            print(f"\nProcessing Latency - Min: {min_processing_latency:.6f}, Max: {max_processing_latency:.6f}")
            print(f"Link Latency - Min: {min_link_latency:.6f}, Max: {max_link_latency:.6f}")
            print(f"Total Latency - Min: {min_total_latency:.6f}, Max: {max_total_latency:.6f}")
            
            # Calculate latency breakdown percentages
            if avg_total_latency > 0:
                processing_percentage = (avg_processing_latency / avg_total_latency) * 100
                link_percentage = (avg_link_latency / avg_total_latency) * 100
                print(f"\nLatency Breakdown:")
                print(f"Processing latency: {processing_percentage:.1f}% of total")
                print(f"Link latency: {link_percentage:.1f}% of total")
            
            print('-------------------------------------')

        return entropy, entropy_mean, entropy_median, entropy_q25
