nextflow.enable.dsl=2

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
        --output_name CEDAR_GRCh37_genotyped\
        --outdir CEDAR

    Mandatory arguments:
      --bfile                       Path to the raw genotypes in plink format (Typically on GRCh37 build)
      --output_name                 Prefix for the output files

    Genotype harmonisation & QC:
      --chain_file                  CrossMap chain file used to convert raw genotype from the source assembly to target assembly of the reference panel (typically GRCh37_to_GRCh38.chain)
      --ref_panel                   VCF file containing the positions, REF and ALT alleles of the imputations reference panel. Used by Genotype Harmonizer (typically 1000 Genomes 30x on GRCh38).
      --ref_genome                  Target reference genome fasta file for CrossMap (typically GRCh38).
      --skip_crossmap               Skip the CrossMap step. Only select this option if the input data uses GRCh38 coordinates (deault: false).

    Phasing & Imputation:
      --eagle_genetic_map           Eagle genetic map file (GRCh38)
      --eagle_phasing_reference     Phasing reference panel for Eagle (typically 1000 Genomes 30x on GRCh38)
      --impute_PAR                  Whether or not to impute the pseudo-autosomal region (PAR) of the X chromosome
      --impute_non_PAR              Whether or not to impute the non-PAR part of the X chromsome
      --minimac_imputation_reference Imputation reference panel for Minimac4 in M3VCF format (typically 1000 Genomes 30x on GRCh38)
      --r2_thresh                   Imputation quality score threshold for filtering poorly imputed variants

    Annotation:
      --annotation_file_23_to_X     Annotation file for bcftools
    
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    
    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
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
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

// Reference panels

Channel
    .from(params.ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .set { ref_panel_ch} 

//Reference genomes
Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

// Eagle
Channel
    .fromPath(params.eagle_genetic_map)
    .ifEmpty { exit 1, "Eagle genetic map file not found: ${params.eagle_genetic_map}" } 
    .set { eagle_genetic_map_ch }

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

// Minimac4
Channel
    .fromPath( "${params.minimac_imputation_reference}*" )
    .ifEmpty { exit 1, "Minimac4 imputation reference not found: ${params.minimac_imputation_reference}" }
    .set { minimac_ref_ch }

// CrossMap
Channel
    .fromPath(params.chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${params.chain_file}" } 
    .set { chain_file_ch }

// bcftools annotation
Channel
    .fromPath(params.annotation_file_23_to_X)
    .ifEmpty { exit 1, "Chromosome annotation file not found: ${params.annotation_file_23_to_X}" }
    .set { annotation_file_23_to_X_ch }



//Chromosome channel
if (params.impute_PAR)
  chromosome_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X")
else
  chromosome_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

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
summary['Harmonisation ref panel']  = params.ref_panel
summary['Eagle genetic map']        = params.eagle_genetic_map
summary['Eagle reference panel']    = params.eagle_phasing_reference
summary['Impute PAR']               = params.impute_PAR
summary['Impute non-PAR']           = params.impute_non_PAR
summary['Minimac4 reference panel'] = params.minimac_imputation_reference
summary['CrossMap chain file']      = params.chain_file
summary['R2 thresh']                = params.r2_thresh
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


include { GenotypeHarmonizer } from './modules/GenotypeHarmonizer'
include { plink_to_vcf } from './modules/preimpute_QC'
include { vcf_fixref; filter_preimpute_vcf; split_by_chr } from './modules/preimpute_QC'
include { annotate } from './modules/preimpute_QC'
include { CrossMap; CrossMap_QC } from './modules/CrossMap'
include { eagle_prephasing } from './modules/eagle'
include { minimac_imputation } from './modules/minimac'
include { filter_vcf; merge_unfiltered_vcf } from './modules/postimpute_QC'
include { initial_filter; preimpute_filter; female_qc; impute_sex; extract_female_samples; split_male_female_dosage; double_male_dosage; merge_male_female} from './modules/impute_non_PAR.nf'

workflow main_flow{
  take:
    bfile
  
  main:
  //Convert input genotypes to GRCh38 coordinates
  if (params.skip_crossmap){
      GenotypeHarmonizer(bfile, ref_panel_ch.collect())
  } else {
      CrossMap(bfile, chain_file_ch.collect())
      GenotypeHarmonizer(CrossMap.out, ref_panel_ch.collect())
  }
  plink_to_vcf(GenotypeHarmonizer.out)
  annotate(plink_to_vcf.out, annotation_file_23_to_X_ch.collect())
  vcf_fixref(annotate.out, ref_genome_ch.collect(), ref_panel_ch.collect())
  filter_preimpute_vcf(vcf_fixref.out)

  //Run imputation on each chromosome
  split_by_chr(filter_preimpute_vcf.out, chromosome_ch)  
  eagle_prephasing(split_by_chr.out, eagle_genetic_map_ch.collect(), phasing_ref_ch.collect())
  minimac_imputation(eagle_prephasing.out, minimac_ref_ch.collect())

  emit:
  minimac_imputation.out[0]
  minimac_imputation.out[1]
  
}

workflow impute_non_PAR{
  take:
    bfile

  main:
  impute_sex(bfile)
  extract_female_samples(impute_sex.out)
  if (params.skip_crossmap){
    GenotypeHarmonizer(impute_sex.out, ref_panel_ch.collect())
  } else {
    CrossMap(impute_sex.out, chain_file_ch.collect())
    GenotypeHarmonizer(CrossMap.out, ref_panel_ch.collect())
  }
  plink_to_vcf(GenotypeHarmonizer.out)
  annotate(plink_to_vcf.out, annotation_file_23_to_X_ch.collect())
  vcf_fixref(annotate.out, ref_genome_ch.collect(), ref_panel_ch.collect())
  initial_filter(vcf_fixref.out)
  female_qc(extract_female_samples.out[0], initial_filter.out)
  preimpute_filter(female_qc.out, initial_filter.out)
  eagle_prephasing(preimpute_filter.out, eagle_genetic_map_ch.collect(), Channel.fromPath(params.eagle_phasing_reference + '/NonPAR/chrX.bcf*').collect())
  minimac_imputation(eagle_prephasing.out, file(params.minimac_imputation_reference + '/NonPAR/chrX.m3vcf.gz'))
  split_male_female_dosage(minimac_imputation.out, extract_female_samples.out[0], extract_female_samples.out[1])
  double_male_dosage(split_male_female_dosage.out[0])
  merge_male_female(double_male_dosage.out, split_male_female_dosage.out[1])

  emit:
  merge_male_female.out[0]
  merge_male_female.out[1]
}

workflow {
  if (params.impute_non_PAR){
    impute_non_PAR(bfile_ch)
    main_flow(bfile_ch)
    merge_unfiltered_vcf(main_flow.out[0].concat(impute_non_PAR.out[0]).collect(), main_flow.out[1].concat(impute_non_PAR.out[1]).collect())
  }
  else{
    main_flow(bfile_ch)
    merge_unfiltered_vcf(main_flow.out[0].collect(), main_flow.out[1].collect())
  }
    filter_vcf(merge_unfiltered_vcf.out)
}