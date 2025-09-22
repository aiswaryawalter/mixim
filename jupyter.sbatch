#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH -J jupyter_gpu
#SBATCH -e jup-%j.err

module load nvidia/latest

IMG=pytorch_25.02-py3.simg

# <-- note the “--” here -->
rz-singularity $IMG -- jupyter-notebook \
     --no-browser --ip="*" --port=8889
