#!/bin/bash

#SBATCH --job-name='chipimputation'
#SBATCH --cpus-per-task=2
#SBATCH --mem=50GB
#SBATCH --output=chipimputation-%j-stdout.log
#SBATCH --error=chipimputation-%j-stderr.log
#SBATCH --time=14-00:00:00


echo "Submitting SLURM job"
nextflow -c /scratch3/users/nanje/chipimputation_evaluate_chips/conf/ilifu/process_ref_h3a_mis_v4.config \
    run /scratch3/users/nanje/chipimputation_evaluate_chips/process_ref_h3a_mis.nf -profile singularity,slurm -resume
