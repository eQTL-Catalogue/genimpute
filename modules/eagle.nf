process eagle_prephasing{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple val(chromosome), file(vcf)
    file eagle_genetic_map
    file phasing_reference

    output:
    tuple val(chromosome), file("chr${chromosome}.phased.vcf.gz")

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=chr${chromosome}.bcf \
    --geneticMapFile=${eagle_genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr${chromosome}.phased \
    --numThreads=8
    """
}

process find_problematic_variants{
    publishDir "${params.outdir}/problem_variants/", mode: 'copy', pattern: "*problem_variants.txt"
    container = 'quay.io/eqtlcatalogue/susier:v21.10.2'

    input:
    file chrX_phased_vcf

    output:
    file("problem_variants.txt")

    script:
    """
    Rscript $baseDir/bin/find_problematic_variants.R --filename ${chrX_phased_vcf} 

    """
}