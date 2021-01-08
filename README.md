# eQTL-Catalogue/genimpute workflow
Genotype imputation and quality control workflow used by the eQTL Catalogue. Performs the following main steps:

**Pre-imputation QC:**
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/). 
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).

**Imputation:**
- Genotype pre-phasing with Egale 2.4.1 
- Genotype imputation with Minimac4

**Post-imputation QC:**
- Convert genotypes to GRCh38 coordinates with CrossMap.py v0.4.1
- Exclude variants with imputation R2 < 0.2
- Keep variants on chromosomes 1-22 and X
- Keep variants with MAF > 0.01

## Running the workflow at UT HPC

```
module load java-1.8.0_40
module load singularity/3.5.3

nextflow run main.nf -profile eqtl_catalogue -resume\
  --bfile /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/PLINK_100718_1018/CEDAR\
  --harmonise_genotypes true\
  --output_name CEDAR_GRCh37_genotyped\
  --outdir results_CEDAR
```

To run the pipeline for the first time, you need to log into the stage1 node to start nextflow there so that it is able to build the nextflow container. Subsequently, you can also use the `srun --pty bash` command to run nextflow as an slurm job interactive job. It is recommended to run srun or nexflow inside a [screen](https://linuxize.com/post/how-to-use-linux-screen/) session.

## Contributors
* Kaur Alasoo
* Liina Anette Pärtel
* Mark-Erik Kodar
