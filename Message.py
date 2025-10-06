class Message:
    def __init__(self, id, type, sender, route, delays, link_delays, pr_target, target_bool):
        self.id = "%d_%d" % (sender.id, id)
        self.type = type  # Dummy or Real packet
        self.sender = sender  # sender object
        self.route = route  # e.g. [S1, mix1, mix2, mix3, R5]
        self.delays = delays  # list of delays at 3 nodes
        self.link_delays = link_delays  # transmission delays between nodes [client->mix1, mix1->mix2, ...]
        self.pr_target = pr_target  # probability of this message being the target message
        self.target_bool = target_bool  # True if this message is a target message
        self.time_left = 0
        self.next_hop_index = 1
    
    def get_total_latency(self):
        """Calculate total latency: sum of processing delays + link delays"""
        processing_latency = sum(self.delays) if self.delays else 0
        link_latency = sum(self.link_delays) if self.link_delays else 0
        total_latency = processing_latency + link_latency
        return total_latency
    
    def get_processing_latency(self):
        """Get total processing delay"""
        return sum(self.delays) if self.delays else 0
    
    def get_link_latency(self):
        """Get total link delay"""
        return sum(self.link_delays) if self.link_delays else 0
