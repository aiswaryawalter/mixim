# Mixim
## A simulation framework for mixnets

This project is a fork of the [original Mixim repository](https://gitlab.esat.kuleuven.be/Iness.BenGuirat/mixim), enhanced with additional features to support advanced topology configurations.

## Features

### Topologies
1. **Stratified** - Traditional layered mix network
2. **Free Route** - Fully connected mesh topology  
3. **Barabasi Albert (BA)** - Scale-free network topology
4. **Cyclic Stratified** - Ring-based layered topology
5. **Grid** - 2D grid mesh topology *(New)*
6. **XRD/Cascade** - Cascade-based routing chains

### Routing Strategies
1. **Source** - Complete route determined at source
2. **Hop-by-hop** - Route determined incrementally at each hop
3. **LARMix** - Latency-aware routing for geographically distributed networks

### Mix Types
1. **Pool** - Threshold-based batching with partial flush
2. **Timed** - Time-based periodic flushing
3. **Poisson** - Individual message processing with exponential delays

## Configuration Guide


### Basic Configuration Structure

Create a `ConfigFile.ini` with the following sections:

```ini
[DEFAULT]
n_clients = 200          # Number of clients
lambda_c = 1            # Message generation rate (1/rate_client)
n_hops = 5              # Number of hops in message path

[TOPOLOGY]
type = stratified       # Topology type
fully_connected = True  # Full connectivity between layers
routing = source        # Routing strategy

[MIXING]
mix_type = pool         # Mix node type
threshold = 100         # Pool size threshold
flush_percent = 1.0     # Percentage of pool to flush
mu = 1.0               # Average processing delay at mixes

[THREATMODEL]
corrupt_mixes = 5       # Number of corrupted mixes
balanced_corruption = True
```
### Stratified Topology
```
[TOPOLOGY]
type = stratified
fully_connected = True
routing = source        # or hopbyhop
n_layers = 5
l_mixes_per_layer = 10

[MIXING]
mix_type = pool
threshold = 50
flush_percent = 0.8
```
### Free Route Topology
```
[TOPOLOGY]
type = free route
fully_connected = True
routing = source        # or hopbyhop or larmix
n_layers = 1
l_mixes_per_layer = 20

[MIXING]
mix_type = pool
threshold = 30
flush_percent = 0.6
```
### Barabasi Albert (BA) Topology
```
[TOPOLOGY]
type = ba topology
routing = source
n_layers = 1
l_mixes_per_layer = 15
m_barabasi_mixes = 3    # Number of edges for new nodes

[MIXING]
mix_type = pool
threshold = 20
flush_percent = 0.6
```
### Cyclic Stratified Topology
```
[TOPOLOGY]
type = cyclic_stratified
routing = source
n_layers = 4
l_mixes_per_layer = 8

[MIXING]
mix_type = pool
threshold = 25
flush_percent = 0.8
```
### Grid Topology
```
[TOPOLOGY]
type = grid
routing = source
grid_width = 5          # Grid width (number of columns)
grid_height = 4         # Grid height (number of rows)
n_layers = 1
l_mixes_per_layer = 20  # Must equal grid_width * grid_height

[MIXING]
mix_type = pool
threshold = 15          # Lower threshold for grid
flush_percent = 0.7

[1]---[2]---[3]
 |     |     |
[4]---[5]---[6]
 |     |     |
[7]---[8]---[9]
```

## Routing Strategies

### Source Routing
```
[TOPOLOGY]
routing = source
```
### Hop-by-Hop Routing
```
[TOPOLOGY]
routing = hopbyhop
```
### LARMix Routing
```
[TOPOLOGY]
routing = larmix
# Works with: stratified, free route

[LARMIX]
tau = 0.1              # Latency sensitivity parameter
latency_bound = 2.0    # Maximum acceptable latency
balancing = True       # Enable load balancing
```
## Mix Type Configurations
### Poisson Mix
```
[MIXING]
mix_type = poisson
mu = 1.0               # Mean processing delay
```
### Pool Mix
```
[MIXING]
mix_type = pool
threshold = 50          # Messages needed to trigger flush
flush_percent = 0.6     # Percentage of pool to flush (0.0-1.0)
```
### Timed Mix
```
[MIXING]
mix_type = time
flush_timeout = 5.0     # Flush interval in time units
```