process plink_to_vcf{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(bed), file(bim), file(fam)

    output:
    file "harmonised.vcf.gz"

    script:
    """
    plink2 --bfile ${bed.simpleName} --recode vcf-iid --out ${bed.simpleName}
    bgzip harmonised.vcf
    """
}

process vcf_fixref{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file input_vcf
    file fasta
    tuple file(ref_panel_vcf), file(ref_panel_vcf_index)

    output:
    file "fixref.vcf.gz"

    script:
    """
    bcftools index ${input_vcf}
    bcftools +fixref ${input_vcf} -- -f ${fasta} -i ${ref_panel_vcf} | \
     bcftools norm --check-ref x -f ${fasta} -Oz -o fixref.vcf.gz
    """
}