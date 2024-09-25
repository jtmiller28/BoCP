#!/bin/bash

#SBATCH --job-name=determine-native-range-array               # Job name
#SBATCH --mail-type=FAIL,ARRAY_TASKS     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=/blue/guralnick/millerjared/BoCP/logs/native-range-array_%A-%a.out                 # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem-per-cpu=5gb                   # Job memory request
#SBATCH --time=00-2:00:00               # Time limit days-hrs:min:sec
#SBATCH --array=1001-4000%20                  # Array Range
#SBATCH --account=guralnick             # We're using Soltis
#SBATCH --qos=guralnick-b                # We're using Soltis resources
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/finished-code/r/003D-determine-native-range.R

## example for running: sbatch /blue/guralnick/millerjared/BoCP/finished-code/bash/003D-determine-native-range-scheduler.sh