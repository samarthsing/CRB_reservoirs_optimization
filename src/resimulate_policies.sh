#!/bin/bash
#SBATCH -D /project/quinnlab/ss9vz/CRB_reservoirs_opt/src/# working directory
#SBATCH -o /project/quinnlab/ss9vz/CRB_reservoirs_opt/job_output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 5
#SBATCH --ntasks-per-node 15
#SBATCH -p parallel        				# Queue name "parallel"
#SBATCH -A quinnlab_paid       					# allocation name
#SBATCH -t 1:00:00       					# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=ss9vz@virginia.edu      # address for email notification
#SBATCH --mail-type=ALL                  	# email at Begin and End of job
module load gcc/9.2.0 openmpi/3.1.6

# Your commands go here
# arguments are <seed> <NFE>
mpic++ -std=c++2a HYSSR_CRB_model_sim.cpp -o run_model
srun ./run_model
