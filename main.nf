Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find study file: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.chromosome, file(row.vcf) ]}
    .set { vcf_ch }

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

vcf_ch.view()
phasing_ref_ch.view()



