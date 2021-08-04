process CrossMap{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(bed), file(bim), file(fam)
    file chain_file

    output:
    tuple file("crossmapped_plink.bed"), file("crossmapped_plink.bim"), file("crossmapped_plink.fam")

    shell: 
    //Converts BIM to BED and converts the BED file via CrossMap. 
    //Finds excluded SNPs and removes them from the original plink file. 
    //Then replaces the BIM with CrossMap's output.
    """
    awk '{print \$1,\$4,\$4+1,\$2,\$5,\$6,\$2 "_" \$5 "_" \$6}' ${bed.simpleName}.bim > crossmap_input.bed
    CrossMap.py bed ${chain_file} crossmap_input.bed crossmap_output.bed
    awk '{print \$7}' crossmap_input.bed | sort > input_ids.txt
    awk '{print \$7}' crossmap_output.bed | sort > output_ids.txt
    comm -23 input_ids.txt output_ids.txt | awk '{split(\$0,a,"_"); print a[1]}' > excluded_ids.txt
    plink2 --bfile ${bed.simpleName} --exclude excluded_ids.txt --make-bed --output-chr MT --out crossmapped_plink
    awk -F'\t' 'BEGIN {OFS=FS} {print \$1,\$4,0,\$2,\$5,\$6}' crossmap_output.bed > crossmapped_plink.bim
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
    bcftools index -t ${vcf.simpleName}_QCd.vcf.gz
    """
}