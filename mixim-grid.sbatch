#!/bin/bash
#SBATCH --job-name=mixim_job
#SBATCH --output=Logs/mixim-%j.out
#SBATCH --error=Logs/mixim-%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G

module load python/3.12
mkdir -p Logs

python3.12 -m venv venv_mixim
source venv_mixim/bin/activate

# Add system monitoring
echo "=== SYSTEM INFO ===" 
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "Available memory: $(free -h)"
echo "Python version: $(python3.12 --version)"
echo "Python path: $(which python3.12)"

which python
python --version


pip install --upgrade pip
pip install -r requirements.txt
# export PYTHONPATH=$HOME/.local/lib/python3.12/site-packages:$PYTHONPATH

echo "Starting job at $(date)"
echo "Available memory: $(free -h)"

echo "=== TESTING IMPORTS ==="
python3.12 -c "from main import main; print('Import successful')" || {
    echo "ERROR: Failed to import main"
    exit 1
}

echo "=== STARTING JOB ==="
python3.12 grid_poolsize_hops_analysis.py 2>&1 | tee Logs/mixim-detailed-${SLURM_JOB_ID}.log

EXIT_CODE=$?
echo "=== JOB COMPLETED with exit code: $EXIT_CODE ==="
echo "End time: $(date)"

exit $EXIT_CODE
echo "Final memory usage: $(free -h)"
