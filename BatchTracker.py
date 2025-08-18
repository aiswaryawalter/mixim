
incoming_batches = {}  # {batch_id: {msg_id: send_time}}
outgoing_batches = {}  # {batch_id: {msg_id: recv_time}}
incoming_outgoing_batch_map = {}  # {client_id: {incoming_batch_id: outgoing_batch_id}}
next_incoming_batch_id = 0
next_outgoing_batch_id = 0