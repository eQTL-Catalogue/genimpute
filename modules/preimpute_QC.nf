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

process filter_preimpute_vcf{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}_preimpute.vcf.gz" else null }
    
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file input_vcf

    output:
    tuple file("filtered.vcf.gz"), file("filtered.vcf.gz.csi")

    script:
    """
    #Add tags
    bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz

    #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
     bcftools filter -e "ALT='.'" |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

     #Index the output file
     bcftools index filtered.vcf.gz
    """
}

process calculate_missingness{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "genotypes.imiss") "${params.output_name}.imiss" else null }
    
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(input_vcf), file(input_vcf_index) 

    output:
    file "genotypes.imiss"

    script:
    """
    vcftools --gzvcf ${input_vcf} --missing-indv --out genotypes
    """
}

process split_by_chr{
    publishDir "${params.outdir}/preimpute/split_chr", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }

    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(input_vcf), file(input_vcf_index)
    each chr

    output:
    tuple val(chr), file("chr_${chr}.vcf.gz")

    script:
    """
    bcftools view -r ${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}