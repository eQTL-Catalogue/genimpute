Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find study file: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.chromosome, file(row.vcf) ]}
    .set { vcf_ch }

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

process eagle_prephasing{
    
    input:
    set chromosome, file(vcf) from vcf_ch
    file genetic_map from genetic_map_ch.collect()
    file phasing_reference from phasing_ref_ch.collect()

    output:
    set chromosome, file("chr_${chromosome}.phased.vcf.gz") into phased_vcf_cf

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
    publishDir "${params.outdir}/vcf/", mode: 'copy', pattern: "*.dose.vcf.gz"
    
    input:
    set chromosome, file(vcf) from phased_vcf_cf
    file imputation_reference from imputation_ref_ch.collect()

    output:
    set chromosome, file("chr_${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    script:
    """
    minimac4 --refHaps ${chromosome}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome
    """
}


