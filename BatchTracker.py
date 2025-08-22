
from collections import defaultdict, Counter

next_incoming_batch_id = 0
next_outgoing_batch_id = 0
incoming_batches = {} 
outgoing_batches = {}  
incoming_outgoing_batch_map = {} 
valids = []  
batch_prob = {}
out_batch_mapping_count = defaultdict(Counter) 
out_msg_mapping_set = {} 
anonymity_set = {}
anonymity_set_size = {}

def compute_batch_permutations(message):
        global valids
        out_batch_id = message.outgoing_batch_id
        out_msg_id = message.outgoing_msg_id
        out_msg_mapping_set[out_msg_id] = set()
        out_msg_time = outgoing_batches[out_batch_id][out_msg_id]
        print(f"==>> OutMsgID: {out_msg_id}")
        print(f"==>> OutMsgTime: {out_msg_time}")
        print(f"==>> Incoming Batches: {incoming_batches}")

        for in_batch_id in incoming_batches:
            len_in = len(incoming_batches[in_batch_id])
            len_out = len(outgoing_batches[out_batch_id])
            if len_in >= len_out:
                print(f"==>> OutBatchMappingCount[{out_batch_id}]: {in_batch_id} Added ")
                out_batch_mapping_count[out_batch_id][in_batch_id] = 0

        for in_batch_id in out_batch_mapping_count[out_batch_id]:
            for inc_msg_id, time_left in incoming_batches[in_batch_id].items():
                print(f"==>> Time Left[{inc_msg_id}]: {time_left}")
                if time_left < out_msg_time:
                    print(f"==>> OutMsgMappingSet[{out_msg_id}]: {inc_msg_id} Added ")
                    out_msg_mapping_set[out_msg_id].add(inc_msg_id)

        print(f"==>> OutMsgMappingSet[{out_msg_id}]: {out_msg_mapping_set[out_msg_id]}")
        if not valids:
            for inc_msg in out_msg_mapping_set[out_msg_id]:
                valids.append({out_batch_id: [inc_msg]})
                i = batchid(inc_msg)
                out_batch_mapping_count[out_batch_id][i] += 1 
            # print(f"==>> Valids in 1st loop: {valids}")
        else:
            # print(f"==>> Valids in 2nd loop: {valids}")
            # print(f"==>> OutMsgMappingSet[{out_msg_id}]: {out_msg_mapping_set[out_msg_id]}")
            temp_valids = []
            for inc_msg in out_msg_mapping_set[out_msg_id]:
                # print(f"==>> IncMsg in the loop: {inc_msg}")
                i = batchid(inc_msg)  
                j = msgid(inc_msg)  
                for x in valids:
                    new_x = x.copy()
                    count = 0
                    msg_list = new_x.get(out_batch_id, [])
                    if msg_list:
                        for v in range(len(msg_list)):
                            # print(f"==>> Valids X[{out_batch_id}]: {msg_list}")
                            # print(f"==>> Valids X[{out_batch_id}][{v}] inside the loop: {msg_list[v]}")
                            b_id = batchid(msg_list[v])
                            m_id = msgid(msg_list[v])
                            if b_id == i and m_id != j:
                                count += 1
                            else:
                                break
                        if count == len(msg_list):
                            new_msg_list = append_msg(msg_list, inc_msg)
                            new_x[out_batch_id] = new_msg_list
                            temp_valids.append(new_x)
                            # print(f"==>> Updating existing {x[out_batch_id]} --> {new_x[out_batch_id]}")
                            out_batch_mapping_count[out_batch_id][i] += 1
                            count = 0
                        else:
                            count = 0
                    else:
                        for batch_msgs in x.values():
                            if batch_msgs:  
                                if batchid(batch_msgs[0]) != i:
                                    count += 1
                        if count == len(x):
                            new_x[out_batch_id] = [inc_msg]
                            temp_valids.append(new_x)
                            # print(f"==>> Adding New Batch to Valids: {new_x[out_batch_id]}")
                            out_batch_mapping_count[out_batch_id][i] += 1
                            count = 0
                        else:
                            count = 0
            if temp_valids:
                valids = temp_valids
                temp_valids = []
        print(f"==>> Number of Valids: {len(valids)}")
        # print(f"==>> Valids: {valids}")
        for x in valids:
            for out_id in x:
                if out_id != out_batch_id:  
                    inc_msgs = x[out_id]
                    if inc_msgs:  
                        inc_id = batchid(inc_msgs[0])
                        out_batch_mapping_count[out_id][inc_id] += 1
        print(f"==>> OutBatchMappingCount: {out_batch_mapping_count}")
        for out_batch in out_batch_mapping_count:
            if out_batch not in batch_prob:
                batch_prob[out_batch] = {}
            non_zero = {}
            for in_batch, count in out_batch_mapping_count[out_batch].items():
                # print(f"==>> OutBatch: {out_batch}, InBatch: {in_batch} Count: {count}")
                prob = count / len(valids) if len(valids) > 0 else 0
                if prob > 0:
                    non_zero[in_batch] = prob
            if non_zero:
                batch_prob[out_batch] = non_zero
                anonymity_set[out_batch] = set(non_zero.keys())
                anonymity_set_size[out_batch] = len(anonymity_set[out_batch])
                print(f"==>> AnonymitySetSize[{out_batch}]: {anonymity_set_size[out_batch]}")
            else:
                if out_batch in batch_prob:
                    del batch_prob[out_batch]
        print(f"==>> BatchProb: {batch_prob}")
        out_batch_mapping_count.clear()
        batch_prob.clear()

def batchid(msg):
    parts = msg.split('_')
    return int(parts[1])

def msgid(msg):
    parts = msg.split('_')
    return int(parts[2])

def append_msg(batch, msg):
    new_batch = batch[:]  # shallow copy
    new_batch.append(msg)
    return new_batch




