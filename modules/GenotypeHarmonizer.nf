process GenotypeHarmonizer{
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