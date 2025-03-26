#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --tasks-per-node=1
#SBATCH --job-name=POKE_all_snw       
#SBATCH --output=POKE_all_snw.out     
#SBATCH --time=48:00:00        
#SBATCH --ntasks=2 
#SBATCH --cpus-per-task=10


srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript POKE_2022_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript POKE_2021_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript POKE_2020_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript POKE_2019_com_snw.R &
wait