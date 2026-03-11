#!/bin/bash

#SBATCH --job-name=array_thin                # Job name
#SBATCH --mail-type=FAIL,ARRAY_TASKS     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=/blue/guralnick/millerjared/BoCP/logs/thin_array_%A-%a.out                 # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10gb                   # Job memory request
#SBATCH --time=00-2:00:00               # Time limit days-hrs:min:sec
#SBATCH --array=1-10%3                  # Array Range
#SBATCH --account=soltis             # We're using Gurlab
#SBATCH --qos=soltis-b                # We're using Gurlab resources
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/finished-code/r/003D-spThin.R

## example for running: sbatch /blue/guralnick/millerjared/BoCP/finished-code/bash/003D-spThin.sh
