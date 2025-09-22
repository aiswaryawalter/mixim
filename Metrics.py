import pandas as pd


class Metrics:
    def __init__(self):
        self.batch_logs = []

    def add_batch_log(self, out_batch_id, true_in_batch_id, anonymity_set_size, anonymity_set, batch_prob, sim_timestamp=None, utc_timestamp=None, window_index=None, n_clients=None, batch_size=None):
        log_entry = {
            "window_index": window_index,
            "out_batch_id": out_batch_id,
            "true_in_batch_id": true_in_batch_id,
            "correct_batch_prob": batch_prob.get(true_in_batch_id, None),
            "correct_batch_is_highest": batch_prob.get(true_in_batch_id, 0) == max(batch_prob.values()) if batch_prob else None,
            "anonymity_set_size": anonymity_set_size,
            "anonymity_set": anonymity_set,
            "n_clients": n_clients,
            "batch_size": batch_size,
            "batch_prob": batch_prob,
            "sim_timestamp": sim_timestamp,
            "utc_timestamp": utc_timestamp,
        }
        self.batch_logs.append(log_entry)

    def save(self, logDir="Logs/", filename_suffix=""):
        filename = f"{logDir}batch_logs{filename_suffix}.csv"
        pd.DataFrame(self.batch_logs).to_csv(filename, index=False)