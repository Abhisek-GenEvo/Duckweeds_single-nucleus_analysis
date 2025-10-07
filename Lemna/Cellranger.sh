#!/bin/bash

#SBATCH -J Lemna_CR            # Job name
#SBATCH -o Lemna_CR.%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p smp                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH -c 20                    # Total number of cores for the single task
#SBATCH --mem 150G                # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 48:00:00              # Run time (hh:mm:ss)

#SBATCH -A m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=achakrab@uni-mainz.de

## Load all necessary modules if needed
module load bio/CellRanger/7.1.0
cellranger mkref --genome=Lemna_Hap1 --fasta=Lemna_Reference_Genome_Complete.fasta --genes=Lemna_Reference_Genome_Complete.gtf --memgb=150 --nthreads=20
cellranger count --id=Lemna_Rep1 --transcriptome=Lemna_Hap1 --fastqs=fastq_files --sample=Rep1_S1,Rep1_S2,Rep1_S3 --include-introns=true --chemistry=auto --create-bam=false --localcores=20 --localmem=150
cellranger count --id=Lemna_Rep2 --transcriptome=Lemna_Hap1 --fastqs=fastq_files --sample=Rep2_S1,Rep2_S2,Rep2_S3 --include-introns=true --chemistry=auto --create-bam=false --localcores=20 --localmem=150
