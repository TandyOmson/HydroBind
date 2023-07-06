#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -l tmpfs=100G
#$ -N hydrobind
#$ -pe smp 36

export KMP_INIT_AT_FORK=FALSE
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OMP_STACKSIZE=4G

module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh

conda activate hydroBind

~/.conda/envs/hydroBind/bin/python main.py

