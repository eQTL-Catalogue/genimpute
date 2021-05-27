process GenotypeHarmonizer_GRCh37{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(study_name_bed), file(study_name_bim), file(study_name_fam)
    tuple file(vcf_file), file(vcf_file_index)

    output:
    tuple file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam")

    script:
    """
    java -Xmx16g -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
}

process GenotypeHarmonizer_GRCh38{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(input_vcf_file), file(input_vcf_file_index)
    tuple file(ref_vcf_file), file(ref_vcf_file_index)

    output:
    tuple file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam")

    script:
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType VCF\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
}