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
    file(chromsome_names) from chr_names_file_ch.collect()

    output:
    tuple val(chr), file("${chr}.renamed.vcf.gz") into build_minimac_ref_ch, build_beagle_ref_ch, create_eagle_bcf_ch

    script:
    """
    bcftools annotate --rename-chrs ${chromsome_names} ${vcf} -Oz -o ${chr}.renamed.vcf.gz
    """
}

process create_m3vcf{
    publishDir "${params.outdir}/m3vcf/", mode: 'copy', pattern: "*.m3vcf.gz"

    container = "quay.io/eqtlcatalogue/minimac3:v2.0.1"

    input:
    tuple val(chr), file(vcf) from build_minimac_ref_ch

    output:
    tuple val(chr), file("${chr}.m3vcf.gz") into m3vcf_ch

    script:
    """
    minimac3 --refHaps ${vcf} --processReference --prefix ${chr}
    """
}

process create_bref3{
    publishDir "${params.outdir}/bref3/", mode: 'copy', pattern: "*.bref3"

    container = "quay.io/eqtlcatalogue/beagle52:21Apr21.304"

    input:
    tuple val(chr), file(vcf) from build_beagle_ref_ch

    output:
    tuple val(chr), file("${chr}.bref3") into bref3_ch

    script:
    """
    zcat ${vcf} | java -jar /usr/bin/bref3.21Apr21.304.jar > ${chr}.bref3
    """
}

process create_eagle_bcf{
    publishDir "${params.outdir}/eagle_ref/", mode: 'copy', pattern: "*.bcf"

    container = "quay.io/biocontainers/bcftools:1.12--h45bccc9_1"

    input:
    tuple val(chr), file(vcf) from create_eagle_bcf_ch

    output:
    tuple val(chr), file("${chr}.bref3") into bcf_ch

    script:
    """
    bcftools view ${vcf} -Ob -o ${chr}.bcf
    """
}