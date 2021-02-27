#!/bin/bash -
#SBATCH -J TraceHeadOn_0_0_100_275                   # Job Name
#SBATCH -o SpEC.stdout                # Output file name
#SBATCH -e SpEC.stderr                # Error file name
#SBATCH -n 24                  # Number of cores
#SBATCH --ntasks-per-node 24        # number of MPI ranks per node
#SBATCH -t 24:0:00   # Run time
#SBATCH -A sxs                # Account name
#SBATCH --no-requeue

# DO NOT MANUALLY EDIT THIS FILE! See `MakeSubmit.py update -h`.
# This is for submitting a batch job on 'wheeler'.
umask 0022
. bin/this_machine.env || echo 'Not using env file.'
set -x
export PATH=$(pwd -P)/bin:$PATH


EvolveGeodesicsWrapper -a="__EmailAddresses__" -f="__TerminationInfoFile__"
module unload python/anaconda2-4.1.1 && module load python/anaconda3-2019.10
python3 --dir /home/mokounko/TracingMerger/ProcessGeodesics.py __Directory__ --trajectories --xturns
