#!/bin/bash

#SBATCH -J Wolffia_CR            # Job name
#SBATCH -o Wolffia_CR.%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH -c 20                    # Total number of cores for the single task
#SBATCH --mem 150G                # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 24:00:00              # Run time (hh:mm:ss)

#SBATCH -A m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=achakrab@uni-mainz.de

## Load all necessary modules if needed
cellranger-arc mkref --config=Wolffia_Haplotype1.config --memgb=150 --nthreads=20
cellranger-arc count --id=Wolffia_Multiome --reference=Wolffia_Haplotype1 --libraries=libraries.csv --localcores=24 --localmem=120
