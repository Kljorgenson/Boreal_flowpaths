#!/bin/bash

#SBATCH --partition=t1standard
#SBATCH --tasks-per-node=1
#SBATCH --job-name=VAUL_all_snw       
#SBATCH --output=VAUL_all_snw.out     
#SBATCH --time=72:00:00        
#SBATCH --ntasks=4 
#SBATCH --cpus-per-task=10


srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript VAUL_2022_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript VAUL_2021_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript VAUL_2020_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript VAUL_2019_com_snw.R &
wait