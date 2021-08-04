process impute_sex{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple file(bed), file(bim), file(fam)
    
    output:
    tuple file("chr23_noHET.bed"), file("chr23_noHET.bim"), file("chr23_noHET.fam")
    
    script:
    """
    plink2 --bfile ${bed.baseName} --impute-sex --make-bed --out sex_imputed
    plink2 --bfile sex_imputed --chr 23 --make-bed --out chr23
    plink2 --bfile chr23 --split-x b37 no-fail --make-bed --out split_x
    plink2 --bfile split_x --chr 23 --make-bed --out chr23_nonPAR
    plink2 --bfile chr23_nonPAR --make-bed --set-hh-missing --output-chr MT --out chr23_noHET
    
    """
}

process extract_female_samples{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'
    
    input:
    tuple file(bed), file(bim), file(fam)
    
    output:
    file("female_samples.txt")
    file("male_samples.txt")
    
    script:
    """
    plink2 --bfile ${bed.baseName} --filter-females --make-bed --out females_only
    plink2 --bfile ${bed.baseName} --filter-males --make-bed --out males_only
    cut -f1 -d' ' females_only.fam > female_samples.txt
    cut -f1 -d' ' males_only.fam > male_samples.txt
    """
}

process initial_filter{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file input_vcf
    
    output:
    file("filtered_X.vcf.gz")
    
    script:
    """
    bcftools filter -e "ALT='.'" ${input_vcf} | bcftools filter -i 'F_MISSING < 0.05' -Oz -o filtered_X.vcf.gz
    """
}

process female_qc{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(txt)
    file(vcf)

    output:
    file("female_passed_regions.txt")

    script:
    """
    bcftools view -S ${txt} ${vcf} -Oz -o females_only.vcf.gz
    bcftools +fill-tags females_only.vcf.gz -Oz -o females_only_tagged.vcf.gz
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' females_only_tagged.vcf.gz |\
        bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        bcftools filter -e "ALT='.'" |\
        bcftools norm -d all |\
        bcftools norm -m+any |\
        bcftools view -m2 -M2 -Oz -o females_only_filtered.vcf.gz
    bcftools query -f "%CHROM\t%POS\n" females_only_filtered.vcf.gz > female_passed_regions.txt
    """
}

process preimpute_filter{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(passed_regions)
    file(vcf)

    output:
    tuple val("X"), file("full_filtered.vcf.gz")

    script:
    """
    bcftools view -T ${passed_regions} ${vcf} -Oz -o full_filtered.vcf.gz
    """
}

process final_filter{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    tuple val(chromosome), file(vcf)

    output:
    file('completely_filtered.vcf.gz')

    script:
    """
    bcftools filter -i 'INFO/R2 > ${params.r2_thresh} & MAF[0]>0.01' ${vcf} | bcftools sort -Oz -o completely_filtered.vcf.gz
    """
}

process split_male_female_dosage{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(input_vcf)
    file(vcf_index)
    file female_samples
    file male_samples

    output:
    file "male_dosage_only.vcf" 
    file "female_dosage_only.vcf.gz" 

    shell:
    """
    #Keep only dosage field
    bcftools annotate -x ^FORMAT/GT,^FORMAT/DS ${input_vcf} -Oz -o dosage_only.vcf.gz
    #Extract male/female samples
    bcftools view -S ${male_samples} dosage_only.vcf.gz -o male_dosage_only.vcf
    bcftools view -S ${female_samples} dosage_only.vcf.gz -Oz -o female_dosage_only.vcf.gz
    """
}

process double_male_dosage{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(vcf) 

    output:
    file "double_male_dosage.vcf" 

    shell:
    //Doubles DS value for every sample. Expects samples to start from the i-th column and DS to be the second value of the columns (delimiter - ":").
    """
    awk -v OFS='\t' '{if (substr(\$1,1,1)!="#") {for (i=10; i<=NF;i++) { n=split(\$i,a,":"); a[2]=a[2]*2;\$i=a[1];for (j=2; j<=n; j++){\$i=\$i ":"a[j]}}}print \$0}' ${vcf} > double_male_dosage.vcf
    """
}

process merge_male_female{
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'

    input:
    file(male_vcf)
    file(female_vcf)

    output:
    file("chr_X_non_PAR.vcf.gz")
    file("chr_X_non_PAR.vcf.gz.csi")

    shell:
    """
    bgzip ${male_vcf}
    bcftools index ${male_vcf}.gz
    bcftools index ${female_vcf}
    bcftools merge ${female_vcf} ${male_vcf}.gz -Oz -o merged.vcf.gz
    bcftools query -l merged.vcf.gz | sort > samples.txt
    bcftools view -S samples.txt merged.vcf.gz -Oz -o chr_X_non_PAR.vcf.gz
    bcftools index chr_X_non_PAR.vcf.gz
    """
}