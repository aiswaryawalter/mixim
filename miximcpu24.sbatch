#!/bin/bash
#SBATCH --job-name=mixim_job
#SBATCH --output=Logs/mixim-%j.out
#SBATCH --error=Logs/mixim-%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=768G

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

python3.12 main.py 2>&1 | tee Logs/mixim-detailed-${SLURM_JOB_ID}.log

echo "=== JOB COMPLETED ==="
echo "End time: $(date)"
echo "Final memory usage: $(free -h)"
