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
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

// Reference panels
Channel
    .from(params.grch37_ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .set { grch37_ref_panel_ch} 

Channel
    .from(params.grch38_ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .set { grch38_ref_panel_ch} 

//Reference genomes
Channel
    .fromPath(params.grch37_ref_genome)
    .ifEmpty { exit 1, "GRCh38 reference genome file not found: ${params.grch37_ref_genome}" } 
    .set { grch37_genome_ch }

Channel
    .fromPath(params.grch38_ref_genome)
    .ifEmpty { exit 1, "GRCh38 reference genome file not found: ${params.grch38_ref_genome}" } 
    .set { grch38_genome_ch }

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

// Beagle
Channel
    .fromPath( "${params.beagle_genetic_map}*" )
    .ifEmpty { exit 1, "Beagle genetic map not found: ${params.eagle_phasing_reference}" }
    .set { beagle_genetic_map_ch }
Channel
    .fromPath( "${params.beagle_imputation_reference}*" )
    .ifEmpty { exit 1, "Beagle imputation reference not found: ${params.beagle_imputation_reference}" }
    .set { beagle_ref_ch }

// CrossMap
Channel
    .fromPath(params.chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${params.chain_file}" } 
    .set { chain_file_ch }

//Chromsome channel
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


include { GenotypeHarmonizer_GRCh37 } from './modules/GenotypeHarmonizer'
include { GenotypeHarmonizer_GRCh38 } from './modules/GenotypeHarmonizer'
include { plink_to_vcf } from './modules/preimpute_QC'
include { plink_to_vcf as plink_to_vcf_grch38 } from './modules/preimpute_QC'
include { vcf_fixref } from './modules/preimpute_QC'
include { vcf_fixref as vcf_fixref_grch38 } from './modules/preimpute_QC'
include { filter_preimpute_vcf } from './modules/preimpute_QC'
include { split_by_chr } from './modules/preimpute_QC'
include { CrossMap; CrossMap_QC } from './modules/CrossMap'
include { eagle_prephasing } from './modules/eagle'
include { beagle_imputation } from './modules/beagle'


workflow{
  //Convert input genotypes to GRCh38 coordinates
  GenotypeHarmonizer_GRCh37(bfile_ch, grch37_ref_panel_ch.collect())
  plink_to_vcf(GenotypeHarmonizer_GRCh37.out)
  vcf_fixref(plink_to_vcf.out, grch37_genome_ch.collect(), grch37_ref_panel_ch.collect())
  CrossMap(vcf_fixref.out, chain_file_ch.collect(), grch38_genome_ch.collect())
  CrossMap_QC(CrossMap.out)
  
  //Perform pre-imputation QC
  GenotypeHarmonizer_GRCh38(CrossMap_QC.out, grch38_ref_panel_ch.collect())
  plink_to_vcf_grch38(GenotypeHarmonizer_GRCh38.out)
  vcf_fixref_grch38(plink_to_vcf_grch38.out, grch38_genome_ch.collect(), grch38_ref_panel_ch.collect())
  filter_preimpute_vcf(vcf_fixref_grch38.out)

  //Run imputation on each chromosome
  split_by_chr(filter_preimpute_vcf.out)
  eagle_prephasing(split_by_chr.out, eagle_genetic_map_ch.collect(), phasing_ref_ch.collect())
  beagle_imputation(eagle_prephasing.out, beagle_genetic_map_ch.collect(), beagle_ref_ch.collect())

}
