process filter_vcf{
  publishDir "${params.outdir}/minimac_out/filtered", mode: 'copy', pattern: "*.vcf.gz*"
  container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(vcf)
    file(index)

    output:
    file("${params.output_name}_filtered.vcf.gz")
    file("${params.output_name}_filtered.vcf.gz.csi")

    shell:
    """
    bcftools filter -i 'INFO/R2 > ${params.r2_thresh}' ${vcf} | \
      bcftools filter -i 'MAF[0] > 0.01' -Oz -o ${params.output_name}_filtered.vcf.gz
    bcftools index ${params.output_name}_filtered.vcf.gz
    """
}

process merge_unfiltered_vcf{
    publishDir "${params.outdir}/minimac_out/unfiltered", mode: 'copy', pattern: "*.vcf.gz*"
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file input_files
    file indices

    output:
    file("${params.output_name}.all_variants.vcf.gz")
    file("${params.output_name}.all_variants.vcf.gz.csi")


    shell:
    """
    bcftools concat ${input_files.join(' ')} -a -Ou | \
    bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' | \
    bcftools sort -Oz -o ${params.output_name}.all_variants.vcf.gz
    bcftools index ${params.output_name}.all_variants.vcf.gz
    """
}