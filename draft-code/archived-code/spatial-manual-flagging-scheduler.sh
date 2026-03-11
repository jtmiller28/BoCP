#!/bin/bash

#SBATCH --job-name=spatial-manual-flagging     # Job name
#SBATCH --mail-type=ALL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=%j.log                  # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=300gb                       # Job memory request
#SBATCH --time=00-48:00:00               # Time limit days-hrs:min:sec
#SBATCH --account=soltis
#SBATCH --qos=soltis-b
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/draft-code/archived-code/manual-spatial-subsetting.R

## example for running: sbatch /blue/guralnick/millerjared/BoCP/draft-code/archived-code/spatial-manual-flagging.sh
