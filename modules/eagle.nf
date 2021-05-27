process eagle_prephasing{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple chromosome, file(vcf)
    file eagle_genetic_map
    file phasing_reference

    output:
    tuple chromosome, file("chr${chromosome}.phased.vcf.gz")

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=ALL.chr${chromosome}.phase3_integrated.20130502.genotypes.bcf \
    --geneticMapFile=${eagle_genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr${chromosome}.phased \
    --numThreads=8
    """
}