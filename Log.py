
class Log:
    def __init__(self):
        self.sent_messages = {
            "MessageID": [], 
            "MessageType": [], 
            "MessageTimeLeft" :[], 
            "MessageDelay": [], 
            "MessageRoute" :[],
            "MessageLinkDelays": [],  # NEW: Store link delays
            "MessageProcessingLatency": [],  # NEW: Sum of processing delays
            "MessageLinkLatency": [],  # NEW: Sum of link delays
            "MessageTotalLatency": [] # NEW: Total latency (processing + link)
            }
        self.received_messages = {
            "MessageID": [], 
            "MessageType": [], 
            "MessageTimeLeft" :[],
            "MessageTimeReceived":[], 
            "MessageDelay": [], 
            "MessageRoute" :[],
            "MessageTarget" : [],
            "MessageLinkDelays": [],  # NEW: Store link delays
            "MessageProcessingLatency": [],  # NEW: Sum of processing delays
            "MessageLinkLatency": [],  # NEW: Sum of link delays
            "MessageTotalLatency": []  # NEW: Total latency (processing + link)
        }
        self.dummy_messages = {
            "DroppingNode":[],
            "DummyID": [], 
            "DummyType": [], 
            "DummyTimeLeft" :[],
            "DummyDelay": [], 
            "DummyRoute" :[], 
            "DummyPr":[],
            "DummyLinkDelays": [],  # NEW: Store link delays
            "DummyProcessingLatency": [],  # NEW: Sum of processing delays
            "DummyLinkLatency": [],  # NEW: Sum of link delays
            "DummyTotalLatency": []  # NEW: Total latency (processing + link)
        }
        self.mix_loads = {
            "MixID": [],
            "MixLayer": [],
            "MessagesProcessed": [],
            "LoadPercentage": [],  # messages processed by this mix / total messages in layer
            "LayerTotalMessages": [],
        }
        self.mix_sending_stats = {
            "MixID": [],
            "MessagesSent": [],
            "MessagesProcessed": [],
            "SendPercentage": [],
            "GridRow": [],
            "GridCol": []
        }
    def dummies_dropped_end_link(self, dummy, dropping_node):
        processing_latency = dummy.get_processing_latency()
        link_latency = dummy.get_link_latency()
        total_latency = dummy.get_total_latency()

        self.dummy_messages["DroppingNode"].append(dropping_node)
        self.dummy_messages["DummyID"].append(dummy.id)
        self.dummy_messages["DummyType"].append(dummy.type)
        self.dummy_messages["DummyTimeLeft"].append(dummy.time_left)
        self.dummy_messages["DummyDelay"].append(dummy.delays)
        self.dummy_messages["DummyRoute"].append(dummy.route)
        self.dummy_messages["DummyPr"].append(dummy.pr_target)
        self.dummy_messages["DummyLinkDelays"].append(dummy.link_delays)
        self.dummy_messages["DummyProcessingLatency"].append(processing_latency)
        self.dummy_messages["DummyLinkLatency"].append(link_latency)
        self.dummy_messages["DummyTotalLatency"].append(total_latency)

    def sent_messages_f(self, msg):
        processing_latency = msg.get_processing_latency()
        link_latency = msg.get_link_latency()
        total_latency = msg.get_total_latency()

        self.sent_messages["MessageID"].append(msg.id)
        self.sent_messages["MessageType"].append(msg.type)
        self.sent_messages["MessageTimeLeft"].append(msg.time_left)
        self.sent_messages["MessageDelay"].append(msg.delays)
        self.sent_messages["MessageRoute"].append(msg.route)
        self.sent_messages["MessageLinkDelays"].append(msg.link_delays)
        self.sent_messages["MessageProcessingLatency"].append(processing_latency)
        self.sent_messages["MessageLinkLatency"].append(link_latency)
        self.sent_messages["MessageTotalLatency"].append(total_latency)

    def received_messages_f(self, msg):
         # Calculate latencies
        processing_latency = msg.get_processing_latency()
        link_latency = msg.get_link_latency()
        total_latency = msg.get_total_latency()

        self.received_messages["MessageID"].append(msg.id)
        self.received_messages["MessageType"].append(msg.type)
        self.received_messages["MessageTimeLeft"].append(msg.time_left)
        self.received_messages["MessageTimeReceived"].append(msg.timeReceived)
        self.received_messages["MessageDelay"].append(msg.delays)
        self.received_messages["MessageRoute"].append(msg.route)
        self.received_messages["MessageTarget"].append(msg.pr_target)
        self.received_messages["MessageLinkDelays"].append(msg.link_delays)
        self.received_messages["MessageProcessingLatency"].append(processing_latency)
        self.received_messages["MessageLinkLatency"].append(link_latency)
        self.received_messages["MessageTotalLatency"].append(total_latency)
    
    def log_mix_load(self, mix_id, layer, messages_processed, layer_total, load_percentage, is_corrupt):
        """Log load information for a mix node"""
        self.mix_loads["MixID"].append(mix_id)
        self.mix_loads["MixLayer"].append(layer)
        self.mix_loads["MessagesProcessed"].append(messages_processed)
        self.mix_loads["LayerTotalMessages"].append(layer_total)
        self.mix_loads["LoadPercentage"].append(load_percentage)

    def log_mix_sending_stats(self, mix_id, messages_sent, messages_processed, send_percentage, grid_row, grid_col):
        """Log mix sending statistics"""
        self.mix_sending_stats["MixID"].append(mix_id)
        self.mix_sending_stats["MessagesSent"].append(messages_sent)
        self.mix_sending_stats["MessagesProcessed"].append(messages_processed)
        self.mix_sending_stats["SendPercentage"].append(send_percentage)
        self.mix_sending_stats["GridRow"].append(grid_row)
        self.mix_sending_stats["GridCol"].append(grid_col)
    
