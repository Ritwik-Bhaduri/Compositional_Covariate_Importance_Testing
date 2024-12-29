#!/bin/bash
#SBATCH -J CRT.eff  # A single job name for the array
#SBATCH -p serial_requeue,shared,janson_cascade,janson,janson_bigmem # Partition
#SBATCH -c 1 # number of cores
#SBATCH -t 0-08:00  # Running time in the format - D-HH:MM
#SBATCH --mem 8000 # Memory request - 1000 corresponds to 1GB
#SBATCH -o ./log_files/output_file_%A_%a.out # Standard output
#SBATCH -e ./log_files/error_file_%A_%a.err # Standard error
Rscript Run.R #SLURM_ARRAY_TASK_ID
#SBATCH --mail-type=END #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ritwikbhaduri@g.harvard.edu #Email to which notifications will be sent

