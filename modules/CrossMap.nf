process CrossMap{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'
    
    input:
    file vcf
    file chain_file
    file target_ref

    output:
    file "${vcf.simpleName}_mapped.vcf.gz"

    shell:
    """
    #Exlcude structural variants, beause they break latest version of CrossMap.py
    bcftools view --exclude-types other ${vcf} -Oz -o ${vcf.simpleName}_noSVs.vcf.gz
    
    #Run CrossMap.py
    CrossMap.py vcf ${chain_file} ${vcf.simpleName}_noSVs.vcf.gz ${target_ref} ${vcf.simpleName}_mapped.vcf
    bgzip ${vcf.simpleName}_mapped.vcf
    """
}

process CrossMap_QC{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file vcf

    output:
    tuple file("${vcf.simpleName}_QCd.vcf.gz"), file("${vcf.simpleName}_QCd.vcf.gz.tbi")

    shell:
    """
    bcftools sort ${vcf} -Oz -o ${vcf.simpleName}_QCd.vcf.gz
    bcftools index ${vcf.simpleName}_QCd.vcf.gz
    """
}