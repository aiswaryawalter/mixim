#!/bin/bash
#SBATCH --job-name=mixim_job
#SBATCH --output=Logs/mixim-%j.out
#SBATCH --error=Logs/mixim-%j.err
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
#SBATCH --mem=16G

module load python/3.12
module load nvidia/latest
mkdir -p Logs

which python
python --version

python3.12 -m pip install --user -r requirements.txt
export PYTHONPATH=$HOME/.local/lib/python3.12/site-packages:$PYTHONPATH

python3.12 main.py
