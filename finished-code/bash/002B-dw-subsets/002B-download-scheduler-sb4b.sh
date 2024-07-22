#!/bin/bash

#SBATCH --job-name=gbif-occurrence-downloads-sb4b      # Job name
#SBATCH --mail-type=ALL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=gbif-occurrence-downloads-sb4b_%j.log                  # Custom output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=10gb                       # Job memory request
#SBATCH --time=00-96:00:00               # Time limit days-hrs:min:sec
#SBATCH --account=guralnick             # We're using Gurlab
#SBATCH --qos=guralnick-b                 # We're using Gurlab resources
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/finished-code/r/002B-dw-subsets/002B-download-occur-data-sb4b.R

## example for running: sbatch /blue/guralnick/millerjared/BoCP/finished-code/bash/002B-dw-subsets/002B-download-scheduler-sb4b.sh
