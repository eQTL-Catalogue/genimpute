//molecular trait data input data
Channel.fromPath(params.vcf_list)
    .ifEmpty { error "Cannot find vcf list file in: ${params.vcf_list}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.chromosome, file(row.path) ]}
    .set { vcf_file_ch }

Channel
    .fromPath(params.chromosome_names)
    .ifEmpty { exit 1, "Chromosome names file not found: ${params.chromosome_names}" } 
    .set { chr_names_file_ch }

process rename_chromosomes{
    container = 'quay.io/biocontainers/bcftools:1.12--h45bccc9_1'

    input:
    tuple val(chr), file(vcf) from vcf_file_ch
    file(chromsome_names) from chr_names_file_ch

    output:
    tuple val(chr), file("${chr}.renamed.vcf.gz") into renamed_vcf_ch

    script:
    """
    bcftools annotate --rename-chrs ${chromsome_names} ${vcf} -Oz -o ${chr}.renamed.vcf.gz
    """
}

process create_m3vcf{
    container = "quay.io/eqtlcatalogue/minimac3:v2.0.1"

    input:
    tuple val(chr), file(vcf) from renamed_vcf_ch

    output:
    tuple val(chr), file("${chr}.m3vcf.gz") into m3vcf_ch

    script:
    """
    minimac3 --refHaps ${vcf} --processReference --prefix ${chr}
    """
}