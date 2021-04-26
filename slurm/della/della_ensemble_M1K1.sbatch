#!/bin/bash
#SBATCH --job-name=M1K1_ens      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=60G        # memory per cpu-core (4G is default)
#SBATCH --time=06:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=mk28@princeton.edu

module purge
module load julia/1.5.0

cd /scratch/gpfs/mk28/CooperationPolarization2/
julia --project run/della_ensemble_M1K1.jl
