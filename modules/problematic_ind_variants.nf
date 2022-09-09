
process find_problematic_variants{
    publishDir "${params.outdir}/problematic_variants/", mode: 'copy', pattern: "problematic_variants.txt"
    publishDir "${params.outdir}/problematic_individuals/", mode: 'copy', pattern: "problematic_individuals.txt"

    container = 'quay.io/eqtlcatalogue/susier:v21.10.2'

    input:
    tuple val(chromosome), file(vcf)

    output:
    path("problematic_variants.txt"), emit: problematic_variants
    path("problematic_individuals.txt"), emit: problematic_individuals

    script:
    """
    Rscript $baseDir/bin/find_problematic_individuals_variants.R --filename ${vcf} 

    """
}

process remove_problematic_variants{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(bed), file(bim), file(fam)
    file problematic_variants

    output:
    tuple file("*_removed_prob_vars.bed"), file("*_removed_prob_vars.bim"), file("*_removed_prob_vars.fam") , emit: bfiles_wo_problematic_variants 

    shell: 
    """
    plink2 --bfile ${bed.simpleName} --exclude ${problematic_variants} --make-bed --output-chr MT --out ${bed.simpleName}_removed_prob_vars
    """
}