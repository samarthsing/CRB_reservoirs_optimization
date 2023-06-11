#!/bin/bash
#SBATCH -D /project/quinnlab/ss9vz/CRB_reservoirs_opt/src/# working directory
#SBATCH -o /project/quinnlab/ss9vz/CRB_reservoirs_opt/job_output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -p dev        				# Queue name "parallel"
#SBATCH -A quinnlab_paid       					# allocation name
#SBATCH -t 1:00:00       					# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=ss9vz@virginia.edu      # address for email notification
#SBATCH --mail-type=ALL                  	# email at Begin and End of job
module load gcc/9.2.0 openmpi/3.1.6

# Your commands go here
# arguments are <seed> <NFE>
g++ -std=c++2a HYSSR_CRB_model_hyssr_copy.cpp -o run_model
srun ./run_model
