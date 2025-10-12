
from collections import defaultdict, Counter
import time, calendar
import logging
import psutil
import sys
import os

# Add logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

next_incoming_batch_id = 0
next_outgoing_batch_id = 0
incoming_batches = {} 
outgoing_batches = {}  
incoming_outgoing_batch_map = {} 
outgoing_to_incoming_batch_map = {}
valids = []  
batch_prob = {}
out_batch_mapping_count = defaultdict(Counter) 
out_msg_mapping_set = {} 
anonymity_set = {}
anonymity_set_size = {}
msg_count = 0
window_size = 1 # number of messages after which to log metrics
window_index = 0
last_metrics_save_time = 0
metrics_save_interval = 0.2 # make it 0.2 later. 1 sim time units; set to 3600 seconds for 1 hour; set to 7200 for 2 hours

def compute_batch_permutations(self, message):
    try:
        batchtracking_start_time = time.time()
        global valids, msg_count, window_index, last_metrics_save_time
        logger.info(f"==>> window_size: {window_size} ===> metrics_save_interval: {metrics_save_interval}")
        msg_count += 1
        out_batch_id = message.outgoing_batch_id
        out_msg_id = message.outgoing_msg_id
        # outgoing_to_incoming_batch_map[out_batch_id] = message.incoming_batch_id
        true_in_batch_id = message.incoming_batch_id
        true_in_msg_id = message.incoming_msg_id
        out_msg_mapping_set[out_msg_id] = set()
        out_msg_time = outgoing_batches[out_batch_id][out_msg_id]
        logger.info(f"==>> OutMsgID: {out_msg_id} ===> IncMsgID: {message.incoming_msg_id}")
        logger.info(f"==>> OutBatchID: {out_batch_id} ===> IncBatchID: {true_in_batch_id}")
        logger.info(f"==>> OutMsgTime: {out_msg_time}")
        logger.info(f"==>> Incoming Batches: {incoming_batches}")

        for in_batch_id in incoming_batches:
            len_in = len(incoming_batches[in_batch_id])
            len_out = len(outgoing_batches[out_batch_id])
            if len_in >= len_out:
                # print(f"==>> OutBatchMappingCount[{out_batch_id}]: {in_batch_id} Added ")
                out_batch_mapping_count[out_batch_id][in_batch_id] = 0

        for in_batch_id in out_batch_mapping_count[out_batch_id]:
            for inc_msg_id, time_left in incoming_batches[in_batch_id].items():
                # print(f"==>> Time Left[{inc_msg_id}]: {time_left}")
                if time_left < out_msg_time:
                    # print(f"==>> OutMsgMappingSet[{out_msg_id}]: {inc_msg_id} Added ")
                    out_msg_mapping_set[out_msg_id].add(inc_msg_id)

        # print(f"==>> OutMsgMappingSet[{out_msg_id}]: {out_msg_mapping_set[out_msg_id]}")
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
                    msg_list = x.get(out_batch_id, [])
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
        logger.info(f"==>> Number of Valids: {len(valids)}")
        # logger.info(f"==>> Valids: {valids}")
        for x in valids:
            for out_id in x:
                if out_id != out_batch_id:  
                    inc_msgs = x[out_id]
                    if inc_msgs:  
                        inc_id = batchid(inc_msgs[0])
                        out_batch_mapping_count[out_id][inc_id] += 1
        logger.info(f"==>> OutBatchMappingCount: {out_batch_mapping_count}")
        for out_batch in out_batch_mapping_count:
            if out_batch not in batch_prob:
                batch_prob[out_batch] = {}
            non_zero = {}
            for in_batch, count in out_batch_mapping_count[out_batch].items():
                # print(f"==>> OutBatch: {out_batch}, InBatch: {in_batch} Count: {count}")
                prob = count / len(valids) if len(valids) > 0 else 0
                if out_batch == out_batch_id and in_batch == true_in_batch_id:
                    logger.info(f"========= Probability of[{out_batch}] of TRUE InBatch [{true_in_batch_id}]: {prob}============")
                if prob > 0:
                    non_zero[in_batch] = prob
            if non_zero:
                batch_prob[out_batch] = non_zero
                anonymity_set[out_batch] = set(non_zero.keys())
                anonymity_set_size[out_batch] = len(anonymity_set[out_batch])
                logger.info(f"==>> AnonymitySetSize[{out_batch}]: {anonymity_set_size[out_batch]}")
            else:
                if out_batch in batch_prob:
                    del batch_prob[out_batch]
        logger.info(f"==>> BatchProb: {batch_prob}")
        if true_in_batch_id not in anonymity_set.get(out_batch_id, set()):
            logger.warning(f"True incoming batch {true_in_batch_id} not in anonymity set for outgoing batch {out_batch_id}")
        # add metrics logging
        utc_timestamp = calendar.timegm(time.gmtime())
        sim_timestamp = self.env.now
        logger.info(f"============ TIME NOW: {sim_timestamp }, UTC (seconds since epoch): {utc_timestamp} ================")
        if msg_count % window_size == 0:
            window_index += 1
            for out_batch in batch_prob:
                self.simulation.Metrics.add_batch_log(
                    out_batch_id=out_batch,
                    true_in_batch_id= outgoing_to_incoming_batch_map.get(out_batch, None),
                    anonymity_set_size=anonymity_set_size.get(out_batch, 0),
                    anonymity_set=anonymity_set.get(out_batch, set()),
                    batch_prob=batch_prob.get(out_batch, {}),
                    sim_timestamp= sim_timestamp, 
                    utc_timestamp=utc_timestamp,
                    window_index=window_index,
                    n_clients=self.simulation.n_clients if hasattr(self.simulation, "n_clients") else None,
                    batch_size=self.simulation.batch_size if hasattr(self.simulation, "batch_size") else None
                )
        # periodic save of metrics
        if sim_timestamp - last_metrics_save_time >= metrics_save_interval:
            job_id = os.environ.get("SLURM_JOB_ID", "")
            filename_suffix = f"_{job_id}_{int(sim_timestamp)}"
            self.simulation.Metrics.save(self.simulation.logDir, filename_suffix)
            last_metrics_save_time = sim_timestamp
        
        # end of metrics logging
        # Batch analysis
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - TotalIncomingBatches: {len(incoming_batches)}")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - TotalOutgoingBatches: {len(outgoing_batches)}")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - OutBatchSize: {len(outgoing_batches.get(out_batch_id, {}))}")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - CandidateInBatches: {len(out_batch_mapping_count[out_batch_id])}")
        batchtracking_duration = time.time() - batchtracking_start_time
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - ProcessingTime: {batchtracking_duration:.4f}s")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]BATCH_ANALYSIS - ValidPermutationsGenerated: {len(valids)}")
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]MEMORY - RSS: {memory_info.rss / 1024 / 1024:.2f} MB")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]MEMORY - VMS: {memory_info.vms / 1024 / 1024:.2f} MB")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]MEMORY - ValidsSizeEstimate: {sys.getsizeof(valids) / 1024:.2f} KB")
        logger.info(f"[{out_msg_id}-{true_in_msg_id}]MEMORY - BatchProbSizeEstimate: {sys.getsizeof(batch_prob) / 1024:.2f} KB")

        # Clear data structures for next message
        out_batch_mapping_count.clear()
        batch_prob.clear()
        anonymity_set.clear()
        anonymity_set_size.clear()
    except Exception as e:
        logger.error(f"Error processing message {out_msg_id}: {str(e)}")
        logger.error(f"Message details: OutBatch={out_batch_id}, InBatch={true_in_batch_id}")
        raise  # Re-raise to not hide the error

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

    




