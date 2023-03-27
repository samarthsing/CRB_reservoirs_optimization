#!/bin/bash
#SBATCH -D /project/quinnlab/ss9vz/CRB_reservoirs_opt/src/
#SBATCH -o /project/quinnlab/ss9vz/CRB_reservoirs_opt/job_output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 22           # Total number of nodes to request (up to 120)
#SBATCH --ntasks-per-node 40           # Number of processors per node (up to 20)
#SBATCH -p parallel           # Queue name "parallel"
#SBATCH -A quinnlab_paid       # allocation name
#SBATCH -t 72:00:00        # Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=ss9vz@virginia.edu              # address for email notification
#SBATCH --mail-type=ALL                  # email at Begin and End of job

#load openmpi module
module load gcc/9.2.0 openmpi/3.1.6

# Your commands go here
# arguments are <seed> <nfe/time> <islands>
srun ./CRB_FCRPS 23 71.2 2
