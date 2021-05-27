process beagle_imputation{
    //publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.dose.vcf.gz"

    container = "quay.io/eqtlcatalogue/beagle52:21Apr21.304"
 
    input:
    tuple val(chromosome), file(vcf)
    file genetic_map
    file imputation_reference

    output:
    tuple val(chromosome), file("chr${chromosome}.imputed.vcf.gz")

    script:
    """
    java -Xss5m -Xmx16g \
      -jar /usr/bin/beagle.21Apr21.304.jar \
        gt=${vcf} \
        ref=chr${chromosome}.bref3 \
        map=plink.chr${chromosome}.GRCh38.map \
        out=chr${chromosome}.imputed \
        chrom=${chromosome}\
        nthreads=8 \
        ne=20000 \
        impute=true \
        gp=true
    """
}