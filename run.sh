#!/bin/bash

#Töö nimi
#SBATCH -J genimpute

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 30:40:00

# mälu
#SBATCH --mem=5G

module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load nextflow
module load squashfs/4.4

nextflow run main.nf \
  -profile tartu_hpc -resume\
  --bfile /gpfs/space/home/a82371/ALSPAC/ALSPAC_ge/ge_merged_plink_files/ALSPAC_merged_w_ge_genos\
  --output_name ALSPAC\
  --outdir ALSPAC\
  --impute_PAR false\
  --chain_file data/NCBI36_to_GRCh38.chain\
  --impute_non_PAR false
