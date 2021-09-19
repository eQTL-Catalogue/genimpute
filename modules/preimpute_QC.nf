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

process annotate {
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file vcf
    file annotation_file

    output:
    file "X_annotated.vcf.gz"

    script:
    """
    bcftools index ${vcf}
    bcftools annotate --rename-chrs ${annotation_file} ${vcf} -Oz -o X_annotated.vcf.gz
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
    tuple file("non_PAR_excluded.vcf.gz"), file("non_PAR_excluded.vcf.gz.csi")

    script:
    """
    #Add tags
    bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz


    #Filter rare and non-HWE variants and those with abnormal alleles and duplicates. 
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
     bcftools filter -e "ALT='.'" |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

    #Index the output file
    bcftools index filtered.vcf.gz


    #Exclude non-PAR region.
    bcftools view filtered.vcf.gz -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X:10001-2781479,X:155701383-156030895 -m2 -M2 -Oz -o non_PAR_excluded.vcf.gz
    bcftools index non_PAR_excluded.vcf.gz
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
    bcftools query -l ${input_vcf} | sort > samples.txt
    bcftools view -r ${chr} -S samples.txt ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}