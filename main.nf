
Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

Channel
    .from(params.ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .into { ref_panel_harmonise_genotypes; ref_panel_vcf_fixref } 

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

Channel
    .fromPath( "${params.minimac_imputation_reference}*" )
    .ifEmpty { exit 1, "Minimac4 imputation reference not found: ${params.minimac_imputation_reference}" }
    .set { imputation_ref_ch }

Channel
    .fromPath(params.eagle_genetic_map)
    .ifEmpty { exit 1, "Eagle genetic map file not found: ${params.eagle_genetic_map}" } 
    .set { genetic_map_ch }

process harmonise_genotypes{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch
    set file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes.collect()

    output:
    set file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam") into harmonised_genotypes

    script:
    if (params.harmonise_genotypes) {
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
    }
    else {
    """
    cp $study_name_bed harmonised.bed
    cp $study_name_bim harmonised.bim
    cp $study_name_fam harmonised.fam
    """
    }
}

process plink_to_vcf{

    input:
    set file(bed), file(bim), file(fam) from harmonised_genotypes

    output:
    file "harmonised.vcf.gz" into harmonised_vcf_ch

    script:
    """
    plink2 --bfile ${bed.simpleName} --recode vcf-iid --out ${bed.simpleName}
    bgzip harmonised.vcf
    """
}

process vcf_fixref{
    
    input:
    file input_vcf from harmonised_vcf_ch
    file fasta from ref_genome_ch.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_vcf_fixref.collect()

    output:
    file "fixref.vcf.gz" into filter_vcf_input

    script:
    """
    bcftools index ${input_vcf}
    bcftools +fixref ${input_vcf} -- -f ${fasta} -i ${vcf_file} | \
     bcftools norm --check-ref x -f ${fasta} -Oz -o fixref.vcf.gz
    """
}

process filter_vcf{

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}.vcf.gz" else null }

    input:
    file input_vcf from filter_vcf_input

    output:
    set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into split_vcf_input, missingness_input

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
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> if (filename == "genotypes.imiss") "${params.output_name}.imiss" else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from missingness_input 

    output:
    file "genotypes.imiss" into missing_individuals

    script:
    """
    vcftools --gzvcf ${input_vcf} --missing-indv --out genotypes
    """
}

process split_by_chr{

    publishDir "genotype_qc/${params.outdir}/by_chr", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from split_vcf_input
    each chr from Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    output:
    tuple val(chr), file("chr_${chr}.vcf.gz") into individual_chromosomes

    script:
    """
    bcftools view -r ${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}

process eagle_prephasing{

    input:
    set chromosome, file(vcf) from individual_chromosomes
    file genetic_map from genetic_map_ch.collect()
    file phasing_reference from phasing_ref_ch.collect()

    output:
    set chromosome, file("chr_${chromosome}.phased.vcf.gz") into phased_vcf_cf

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=ALL.chr${chromosome}.phase3_integrated.20130502.genotypes.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr_${chromosome}.phased \
    --numThreads=8
    """
}

process minimac_imputation{
    publishDir "${params.outdir}/vcf/", mode: 'copy', pattern: "*.dose.vcf.gz"
 
    input:
    set chromosome, file(vcf) from phased_vcf_cf
    file imputation_reference from imputation_ref_ch.collect()

    output:
    set chromosome, file("chr_${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    script:
    """
    minimac4 --refHaps ${chromosome}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome
    """
}


