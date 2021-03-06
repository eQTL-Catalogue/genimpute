def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     eQTL-Catalogue/genimpute v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile eqtl_catalogue -resume\
        --bfile /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/PLINK_100718_1018/CEDAR\
        --harmonise_genotypes true\
        --output_name CEDAR_GRCh37_genotyped\
        --outdir CEDAR

    Mandatory arguments:
      --bfile                       Path to the TSV file containing pipeline inputs (VCF, expression matrix, metadata)
      --output_name                 Prefix for the output files

    Genotype harmonisation & QC:
      --harmonise_genotypes         Run GenotypeHarmonizer on the raw genotypes to correct flipped/swapped alleles (default: true)
      --ref_panel                   Reference panel used by GenotypeHarmonizer. Ideally should match the reference panel used for imputation.
      --ref_genome                  Reference genome fasta file for the raw genotypes (typically GRCh37).

    Phasing & Imputation:
      --eagle_genetic_map           Eagle genetic map file
      --eagle_phasing_reference     Phasing reference panel for Eagle (typically 1000 Genomes Phase 3)
      --minimac_imputation_reference Imputation reference panel for Minimac4 in M3VCF format (typically 1000 Genomes Phase 3)
      --r2_thresh                   Imputation quality score threshold for filtering poorly imputed variants

    CrossMap.py:
      --target_ref                  Reference genome fasta file for the target genome assembly (e.g. GRCh38)
      --chain_file                  Chain file to translate genomic cooridnates from the source assembly to target assembly
    
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    
    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Define input channels
Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

Channel
    .from(params.ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .into { ref_panel_harmonise_genotypes; ref_panel_vcf_fixref } 

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

Channel
    .fromPath( "${params.minimac_imputation_reference}*" )
    .ifEmpty { exit 1, "Minimac4 imputation reference not found: ${params.minimac_imputation_reference}" }
    .set { imputation_ref_ch }

Channel
    .fromPath(params.eagle_genetic_map)
    .ifEmpty { exit 1, "Eagle genetic map file not found: ${params.eagle_genetic_map}" } 
    .set { genetic_map_ch }

Channel
    .fromPath(params.chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${params.chain_file}" } 
    .set { chain_file_ch }

Channel
    .fromPath(params.target_ref)
    .ifEmpty { exit 1, "CrossMap.py target reference genome file: ${params.target_ref}" } 
    .set { target_ref_ch }

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
eQTL-Catalogue/genimpute v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'eQTL-Catalogue/genimpute'
summary['Pipeline Version']         = workflow.manifest.version
summary['Run Name']                 = custom_runName ?: workflow.runName
summary['PLINK bfile']              = params.bfile
summary['Reference genome']         = params.ref_genome
summary['Harmonise genotypes']      = params.harmonise_genotypes
summary['Harmonisation ref panel']  = params.ref_panel
summary['Eagle genetic map']        = params.eagle_genetic_map
summary['Eagle reference panel']    = params.eagle_phasing_reference
summary['Minimac4 reference panel'] = params.minimac_imputation_reference
summary['CrossMap reference genome'] = params.target_ref
summary['CrossMap chain file']      = params.chain_file
summary['R2 thresh']      = params.chain_file
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['Output name']              = params.output_name
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']            = params.awsregion
   summary['AWS Queue']             = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}
process impute_sex{
    input:
    tuple file(bed), file(bim), file(fam) from bfile_ch
    
    output:
    tuple file("chr23_noHET.bed"), file("chr23_noHET.bim"), file("chr23_noHET.fam") into extract_female_ch, harmonize_input_ch
    
    script:
    """
    plink2 --bfile ${bed.baseName} --impute-sex --make-bed --out sex_imputed
    plink2 --bfile sex_imputed --chr 23 --make-bed --out chr23
    plink2 --bfile chr23 --split-x b37 no-fail --make-bed --out split_x
    plink2 --bfile split_x --chr 23 --make-bed --out chr23_noPAR
    plink2 --bfile chr23_noPAR --make-bed --set-hh-missing --out chr23_noHET
    
    """
}

process extract_female_samples{
    
    input:
    tuple file(bed), file(bim), file(fam) from extract_female_ch
    
    output:
    file("female_samples.txt") into female_sample_list_ch
    
    script:
    """
    plink2 --bfile ${bed.baseName} --filter-females --make-bed --out females_only
    cut -f1 -d' ' females_only.fam > female_samples.txt
    """
 }

 process genotype_harmonizer{
    input:
    tuple file(bed), file(bim), file(fam) from harmonize_input_ch
    tuple file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes.collect()

    output:
    tuple file("harmonized.bed"), file("harmonized.bim"), file("harmonized.fam") into harmonised_genotypes
   
    script:
    """
    sed 's/^23/X/g' ${bim.baseName}.bim > new.bim
    mv new.bim ${bim.baseName}.bim
    
    module load java-1.8.0_40
    java -jar ~/software/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar\
     --input CEDAR_chr23_noHET\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonized
    """
 }
 
 process process make_vcf{
    input:
    tuple file(bed), file(bim), file(fam) from harmonised_genotypes
    file fasta from ref_genome_ch.collect()
    tuple file(vcf_file), file(vcf_file_index) from ref_panel_vcf_fixref.collect()
    
    output:
    file("harmonised_chrX_fixref.vcf.gz") into filter_vcf_input
    
    script:
    """
    plink2 --bfile ${bed.baseName} --recode vcf-iid --out harmonised_vcf
    printf '23\\tX\\n' > 23_to_x.tsv
    bcftools annotate --rename-chrs 23_to_x.tsv harmonised_vcf.vcf -Oz -o harmonised_chrX.vcf.gz
    bcftools filter -e "ALT='.'" harmonised_chrX.vcf.gz | bcftools filter -i 'F_MISSING < 0.05' -Oz -o harmonised_chrX_filtered.vcf.gz
    bcftools +fixref harmonised_chrX_filtered.vcf.gz -Oz -o harmonised_chrX_fixref.vcf.gz -- \
    -f ${fasta} \
    -i ${vcf_file}
    """
  }

process female_qc{
    input:
    file(txt) from female_sample_list_ch
    file(vcf) from filter_vcf_input

    output:
    file(female_passed_regions.txt) into passed_regions

    script:
    """
    bcftools view -S ${txt} ${vcf} -Oz -o females_only.vcf.gz
    bcftools +fill-tags females_only.vcf.gz -Oz -o females_only_tagged.vcf.gz
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' females_only_tagged.vcf.gz |\
        bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        bcftools filter -e "ALT='.'" |\
        bcftools norm -d all |\
        bcftools norm -m+any |\
        bcftools view -m2 -M2 -Oz -o females_only_filtered.vcf.gz
    bcftools query -f "%CHROM\t%POS\n" females_only_filtered.vcf.gz > female_passed_regions.txt
    
    """
}