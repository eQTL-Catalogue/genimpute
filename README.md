# imputation-pipeline
Genotype imputation pipeline using Egale 2.4.1 and Minimac4

## Running the workflow at UT HPC

```
module load java-1.8.0_40
module load singularity/3.5.3

nextflow run main.nf -profile eqtl_catalogue -resume\
  --bfile /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/PLINK_100718_1018/CEDAR\
  --harmonise_genotypes true\
  --output_name CEDAR_GRCh37_genotyped\
  --outdir CEDAR
```

To run the pipeline for the first time, you need to log into the stage1 node to start nextflow there so that it is able to build the nextflow container. Subsequently, you can also use the `srun --pty bash` command to run nextflow as an slurm job interactive job. It is recommended to run srun or nexflow inside a [screen](https://linuxize.com/post/how-to-use-linux-screen/) session.

## Contributors
* Kaur Alasoo
* Liina Anette Pärtel
* Mark-Erik Kodar
