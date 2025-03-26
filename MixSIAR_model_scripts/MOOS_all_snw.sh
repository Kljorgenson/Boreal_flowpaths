#!/bin/bash

#SBATCH --partition=t1standard
#SBATCH --tasks-per-node=1
#SBATCH --job-name=MOOS_all_snw       
#SBATCH --output=MOOS_all_snw.out     
#SBATCH --time=72:00:00        
#SBATCH --ntasks=6 
#SBATCH --cpus-per-task=10


srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2022_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2021_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2020_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2019_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2018_com_snw.R &
srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript MOOS_2015_com_snw.R &
wait