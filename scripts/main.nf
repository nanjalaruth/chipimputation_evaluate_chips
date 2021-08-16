#!/usr/bin/env nextflow

/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Create a test data from a vcf file
 * How to run: nextflow run main.nf --source file.vcf
 *
*/

// Show help emssage
//if (params.help){
//    helpMessage()
//    exit 0
//}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '19.04.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

//source  = Channel.from(params.source)
source  = Channel.fromPath(params.source_dataset).view()
subset_size = params.subset_size

process get_snp_list {
    tag "get_snp_list_${base}_${subset_size}"
//    publishDir "${params.outDir}/testdata", overwrite: true, mode:'copy', pattern: "*map"
    input:
        file(vcf_file) from source
    output:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target) into get_snp_list
        file(chromFile) into chromosome
        file(vcf_file) into vcfFile
    script:
        base = file(vcf_file.baseName).baseName
        chromFile = "${base}.chrom"
        subset_map_ref = "refPanel_testdata_${subset_size}.map"
        subset_map_target = "target_testdata_${subset_size/4}.map"
        """
        bcftools query -f '%CHROM\\t%POS\\t%POS\\n' ${vcf_file} > ${base}.map
        awk -F' ' '{print \$1}' ${base}.map | sort -n | uniq > ${chromFile}
        head -n${subset_size/2} ${base}.map > ${subset_map_ref}
        tail -n ${subset_size/2} ${base}.map >> ${subset_map_ref}
        sort -R ${subset_map_ref} | tail -n ${subset_size/4} > ${subset_map_target}
        """
}

chromosomes = file(chromosome.toSortedList().val[0]).readLines().unique().collect { it as int }.sort()

process subset_vcf {
    tag "subset_vcf_${base}_${subset_size}"
    publishDir "${params.outDir}/testdata", overwrite: true, mode:'copy', pattern: "*.vcf.gz"
    input:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target) from get_snp_list
    output:
        set file(vcf_file), file(subset_map_ref), file(subset_map_target), file(vcf_ref), file(vcf_target) into subset_vcf_ref
    script:
        base = file(vcf_file.baseName).baseName
        vcf_ref = "refPanel_testdata.vcf.gz"
        vcf_target = "target_testdata.vcf.gz"
        """
        tabix ${vcf_file}
        bcftools view --regions-file ${subset_map_ref} ${vcf_file} -Oz -o ${vcf_ref}
        tabix ${vcf_ref}
        bcftools view --regions-file ${subset_map_target} ${vcf_ref} -Oz -o ${vcf_target}
        """
}


process split_vcf_to_chrm {
    tag "split_${base}_${chrm}"
    label "medium"
    input:
        each chrm from chromosomes
        set file(vcf_file), file(subset_map_ref), file(subset_map_target), file(vcf_ref), file(vcf_target) from subset_vcf_ref
    output:
        set chrm, file(vcf_ref_chrm) into split_vcf_to_chrm
    script:
        base = file(vcf_ref.baseName).baseName
        vcf_ref_chrm = "${base}_${chrm}.vcf.gz"
        """
        tabix ${vcf_ref}
        bcftools view \
            --regions ${chrm} \
            -m2 -M2 -v snps \
            ${vcf_ref} \
            -Oz -o ${vcf_ref_chrm}
        """
}


process phase_vcf_chrm {
    tag "phase_${base}_${chrm}"
    input:
        set chrm, file(vcf_file_chrm) from split_vcf_to_chrm
    output:
        set chrm, file(vcf_file_chrm), file("${file_out}.vcf.gz") into phase_vcf_chrm
    script:
        base = file(vcf_file_chrm.baseName).baseName
        file_out = "${base}_phased"
        """
        nblines=\$(zcat ${vcf_file_chrm} | grep -v '^#' | wc -l)
        if (( \$nblines > 0 ))
        then
            tabix ${vcf_file_chrm}
            eagle \
                --vcfTarget=${vcf_file_chrm} \
                --geneticMapFile=${params.eagle_genetic_map} \
                --vcfRef=${params.reference_vcf} \
                --vcfOutFormat=z \
                --chrom=${chrm} \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
            if [ ! -f "${file_out}.vcf.gz" ]; then
                touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
            fi
        else
            touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
        fi
        """
}


process vcf_to_m3vcf {
    tag "m3vcf_${ref}_${chrm}"
    label "medmem"
    publishDir "${params.outDir}/testdata", overwrite: true, mode:'copy'
    input:
        set chrm, file(vcf_chrm), file(vcf_chrm_phased) from phase_vcf_chrm
    output:
        set chrm, file(vcf_chrm), file(vcf_chrm_phased), file(m3vcf_chrm) into vcf_to_m3vcf
    script:
        base = file(vcf_chrm_phased.baseName).baseName
        m3vcf_chrm = "${base}.m3vcf.gz"
        
        """
        Minimac3 \
            --refHaps ${vcf_chrm_phased} \
            --processReference \
            --prefix ${base}
        """
}