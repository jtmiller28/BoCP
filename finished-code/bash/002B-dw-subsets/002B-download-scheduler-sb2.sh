#!/bin/bash

#SBATCH --job-name=gbif-occurrence-downloads-sb2      # Job name
#SBATCH --mail-type=ALL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=%j.log                  # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --constraint=rome                # Constrains to specific node (trying to throw things around the cluster to avoid IP overping)
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=100gb                       # Job memory request
#SBATCH --time=00-96:00:00               # Time limit days-hrs:min:sec
#SBATCH --qos=soltis-b
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/finished-code/r/002B-dw-subsets/002B-download-occur-data-sb2.R

## example for running: sbatch /blue/guralnick/millerjared/BoCP/finished-code/bash/002B-dw-subsets/002B-download-scheduler-sb2.sh
