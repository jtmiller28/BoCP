#!/bin/bash

#SBATCH --job-name=coordinateClean-bigmem             # Job name
#SBATCH --mail-type=FAIL,ARRAY_TASKS     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=/blue/guralnick/millerjared/BoCP/logs/coordinateClean/coordinateCleaning-bigmem%A-%a.out                 # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem-per-cpu=80gb                   # Job memory request
#SBATCH --time=00-96:00:00               # Time limit days-hrs:min:sec
#SBATCH --account=guralnick             # We're using Guralnick
#SBATCH --qos=guralnick-b                # We're using Guralnick resources
pwd; hostname; date

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/BoCP/finished-code/r/007-coordinateClean-bigmem.R

done

date
## example for running: sbatch /blue/guralnick/millerjared/BoCP/finished-code/bash/004-build-species-tables-scheduler.sh