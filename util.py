import numpy as np
import random
import csv

def XRD_New(ListCascades):
    groups_list = []
    group1 = []
    group1.append(ListCascades[1])
    group1.append(ListCascades[2])
    group1.append(ListCascades[3])
    groups_list.append(group1)

    group2 = []
    group2.append(ListCascades[1])
    group2.append(ListCascades[4])
    group2.append(ListCascades[5])
    groups_list.append(group2)

    group3 = []
    group3.append(ListCascades[2])
    group3.append(ListCascades[4])
    group3.append(ListCascades[6])
    groups_list.append(group3)

    group4 = []
    group4.append(ListCascades[3])
    group4.append(ListCascades[5])
    group4.append(ListCascades[6])
    groups_list.append(group4)
    return groups_list

def Weights (Number_layer, Number_Mix_Per_Layer):
    probability = []
    for i in range(Number_layer):
        probabTem = []
        for j in range(Number_Mix_Per_Layer):
            probabTem.append((1) / (Number_Mix_Per_Layer))
        probability.append(probabTem)
    return probability

def load_latency_matrix(filename):
    DEFAULT_INTERSERVER_LAT = 0.05
    latency_matrix = {}
    same_server_count = 0 
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            source_id, target_id, latency = row
            # --- Handle same-server latency ---
            if source_id == target_id:
                latency = DEFAULT_INTERSERVER_LAT
                same_server_count += 1
            # print(f"source_id -> {source_id}  target_id -> {target_id}  latency -> {latency}")
            latency_matrix[(int(source_id), int(target_id))] = float(latency)/1000
        # --- Logging summary ---
        print(f"[LatencyMatrix] Loaded {len(latency_matrix)} entries.")
        print(f"[LatencyMatrix] Same-server pairs adjusted: {same_server_count}")
    return latency_matrix

def load_server_coords(filename):
    node_coords = {}
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            node_id = row['id']
            latitude = float(row['latitude'])
            longitude = float(row['longitude'])
            node_coords[node_id] = (latitude, longitude)
    return node_coords






