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

process harmonise_genotypes{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch
    set file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes.collect()

    output:
    set file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam") into harmonised_genotypes

    script:
    if (params.harmonise_genotypes) {
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
    }
    else {
    """
    cp $study_name_bed harmonised.bed
    cp $study_name_bim harmonised.bim
    cp $study_name_fam harmonised.fam
    """
    }
}

process plink_to_vcf{
    input:
    set file(bed), file(bim), file(fam) from harmonised_genotypes

    output:
    file "harmonised.vcf.gz" into harmonised_vcf_ch

    script:
    """
    plink2 --bfile ${bed.simpleName} --recode vcf-iid --out ${bed.simpleName}
    bgzip harmonised.vcf
    """
}

process vcf_fixref{
    input:
    file input_vcf from harmonised_vcf_ch
    file fasta from ref_genome_ch.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_vcf_fixref.collect()

    output:
    file "fixref.vcf.gz" into filter_vcf_input

    script:
    """
    bcftools index ${input_vcf}
    bcftools +fixref ${input_vcf} -- -f ${fasta} -i ${vcf_file} | \
     bcftools norm --check-ref x -f ${fasta} -Oz -o fixref.vcf.gz
    """
}

process filter_preimpute_vcf{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}_preimpute.vcf.gz" else null }

    input:
    file input_vcf from filter_vcf_input

    output:
    set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into split_vcf_input, missingness_input

    script:
    """
    #Add tags
    bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz

    #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
     bcftools filter -e "ALT='.'" |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

     #Index the output file
     bcftools index filtered.vcf.gz
    """
}

process calculate_missingness{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "genotypes.imiss") "${params.output_name}.imiss" else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from missingness_input 

    output:
    file "genotypes.imiss" into missing_individuals

    script:
    """
    vcftools --gzvcf ${input_vcf} --missing-indv --out genotypes
    """
}

process split_by_chr{
    publishDir "${params.outdir}/preimpute/split_chr", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }
    
    input:
    tuple file(input_vcf), file(input_vcf_index) from split_vcf_input
    each chr from Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    output:
    tuple val(chr), file("chr_${chr}.vcf.gz") into individual_chromosomes

    script:
    """
    bcftools view -r ${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}

process eagle_prephasing{
    input:
    tuple chromosome, file(vcf) from individual_chromosomes
    file genetic_map from genetic_map_ch.collect()
    file phasing_reference from phasing_ref_ch.collect()

    output:
    tuple chromosome, file("chr_${chromosome}.phased.vcf.gz") into phased_vcf_cf

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=ALL.chr${chromosome}.phase3_integrated.20130502.genotypes.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr_${chromosome}.phased \
    --numThreads=8
    """
}

process minimac_imputation{
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.dose.vcf.gz"
 
    input:
    set chromosome, file(vcf) from phased_vcf_cf
    file imputation_reference from imputation_ref_ch.collect()

    output:
    tuple chromosome, file("chr_${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    script:
    """
    minimac4 --refHaps ${chromosome}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome
    """
}

process crossmap_genotypes{
    input:
    file chain_file from chain_file_ch.collect()
    file target_ref from target_ref_ch.collect()
    tuple chromosome, file(vcf) from imputed_vcf_cf

    output:
    file "${vcf.simpleName}_mapped.vcf.gz" into mapped_vcf_ch

    shell:
    """
    #Exlcude structural variants, beause they break latest version of CrossMap.py
    bcftools view --exclude-types other ${vcf} -Oz -o ${vcf.simpleName}_noSVs.vcf.gz
    
    #Run CrossMap.py
    CrossMap.py vcf ${chain_file} ${vcf.simpleName}_noSVs.vcf.gz ${target_ref} ${vcf.simpleName}_mapped.vcf
    bgzip ${vcf.simpleName}_mapped.vcf
    """
}

process filter_vcf{
    input:
    each file(vcf) from mapped_vcf_ch
    val r2_thresh from Channel.from(params.r2_thresh)

    output:
    file "${vcf.simpleName}_filtered.vcf.gz" into filtered_vcfs

    shell:
    """
    bcftools filter -i 'INFO/R2 > ${r2_thresh}' ${vcf} -Oz -o ${vcf.simpleName}_filtered.vcf.gz
    """
}

process merge_vcf{
    input:
    file input_files from filtered_vcfs.collect()

    output:
    file "output.vcf.gz" into merged_vcf_ch

    shell:
    """
    bcftools concat ${input_files.join(' ')} | bcftools sort -Oz -o output.vcf.gz
    """
}

process keep_chromosomes{
    publishDir "${params.outdir}/postimpute/crossmap_vcf", mode: 'copy',
        saveAs: {filename -> if (filename == "renamed.vcf.gz") "${params.output_name}.vcf.gz" else null }
    
    input:
    file input_vcf from merged_vcf_ch

    output:
    file "renamed.vcf.gz" into final_vcf_ch

    shell:
    """
    bcftools index ${input_vcf}
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X ${input_vcf} |\
     bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Oz -o renamed.vcf.gz
    """
}

process maf_filter{
    publishDir "${params.outdir}/postimpute/crossmap_vcf", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}.MAF001.vcf.gz" else null }
    
    input:
    file input_vcf from final_vcf_ch

    output:
    file "filtered.vcf.gz" into filtered_vcf_channel

    shell:
    """
    bcftools filter -i 'MAF[0] > 0.01' ${input_vcf} -Oz -o filtered.vcf.gz
    """
}
