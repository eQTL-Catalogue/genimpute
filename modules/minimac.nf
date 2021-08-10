process minimac_imputation{
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.dose.vcf.gz"
 
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple val(chromosome), file(vcf)
    file imputation_reference

    output:
    file("chr_${chromosome}.dose.vcf.gz")
    file("chr_${chromosome}.dose.vcf.gz.csi")

    script:
    """
    minimac4 --refHaps chr${chromosome}.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome \
    --minRatio 0.05
    bcftools index chr_${chromosome}.dose.vcf.gz
    """
}
