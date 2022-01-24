#!/bin/bash

#SBATCH --job-name='chipimputation'
#SBATCH --cpus-per-task=2
#SBATCH --mem=50GB
#SBATCH --output=chipimputation-%j-stdout.log
#SBATCH --error=chipimputation-%j-stderr.log
#SBATCH --time=14-00:00:00


echo "Submitting SLURM job"
nextflow -c process_ref_h3a_mis_v4.config run process_ref_h3a_mis.nf -profile singularity,slurm -resume
