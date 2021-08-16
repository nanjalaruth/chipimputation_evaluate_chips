#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include { filter_min_ac; split_multi_allelic; dataset_qc_dupl; get_chromosome; get_chromosome_vcf; check_chromosome; check_files; check_chromosome_vcf; check_mismatch; no_mismatch ; target_qc as target_qc; target_qc as target_qc1; qc_site_missingness as qc_site_missingness1; qc_site_missingness as qc_site_missingness2; sites_only ; combine_vcfs ; combine_infos; combine_csvs as combine_freqs; combine_vcfs_chrm; } from './modules/qc'


// def check_files(file_list) {
//     file_list.each { myfile ->
//         if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
//     }
// }

process phasing_vcf {
    tag "phase_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf), file(dataset_sample), refpanel, file(refpanel_vcf), file(eagle_genetic_map)

    output:
        tuple chrm, dataset, file("${file_out}.vcf.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        file_out = "${base}_phased"
        """
        tabix ${dataset_vcf}
        tabix ${refpanel_vcf}
        eagle \
            --numThreads=${task.cpus} \
            --vcfTarget=${dataset_vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfRef=${refpanel_vcf} \
            --vcfOutFormat=z \
            --chrom=${chrm} \
            --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
        """
}


process phasing_vcf_no_ref {
    tag "phase_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased", overwrite: true, mode:'copy'
    label "bigmem10"

    input:
        tuple dataset, chrm, file(dataset_vcf), file(eagle_genetic_map)

    output:
        tuple chrm, dataset, file("${file_out}.vcf.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        file_out = "${base}_phased"
        """
        eagle \
            --numThreads=${task.cpus} \
            --vcf=${dataset_vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfOutFormat=z \
            --chrom=${chrm} \
            --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
        tabix "${file_out}.vcf.gz"
        """
}

process tabix_phasing_vcf_no_ref {
    tag "tabix_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased_no_ref", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf)

    output:
        tuple chrm, dataset, file(dataset_vcf), file("${dataset_vcf}.tbi")

    script:
        """
        tabix ${dataset_vcf}
        """
}


process vcf_to_m3vcf {
    tag "m3vcf_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/m3vcfs", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf)

    output:
        tuple chrm, dataset, file(dataset_m3vcf)

    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_m3vcf = "${base}.m3vcf.gz"
        """
        minimac3 \
            --refHaps ${dataset_vcf} \
            --processReference \
            --prefix ${base}
        """
}

process vcf_legend {
    tag "legend_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/legends", mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf)

    output:
        tuple chrm, dataset, file("${base}.legend.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        """
        echo -e 'id position a0 a1 afr.aaf super_pop.aaf afr.maf super_pop.maf' > ${base}.legend
        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' ${dataset_vcf} -Ob -o ${base}_temp1.bcf
        bcftools +fill-tags ${base}_temp1.bcf -Ob -o ${base}_temp2.bcf
        bcftools query -f '%ID %POS %REF %ALT %INFO/AF %INFO/AF %INFO/MAF %INFO/MAF\\n' ${base}_temp2.bcf >> ${base}.legend
        bgzip ${base}.legend
        rm -f ${base}_temp*
        """
}

process vcf_to_bcf {
    tag "bcf_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/bcfs", mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf)

    output:
        tuple chrm, dataset, file("${base}.bcf"), file("${base}.bcf.csi")

    script:
        base = file(dataset_vcf.baseName).baseName
        """
        bcftools view ${dataset_vcf} -Ob -o${base}.bcf
        tabix -f ${base}.bcf
        """
}


workflow{
    check_files([params.eagle_genetic_map])
    params.datasets.each { dataset, dataset_vcf, dataset_sample ->
        datasets_ch = Channel.of(dataset)
            .combine(Channel.fromPath(dataset_vcf))
            .combine(Channel.fromPath(dataset_sample))
    }

    // Get chromosome
    get_chromosome(datasets_ch.map{dataset, vcf, sample -> [dataset, vcf]})

    // Checkk REF mismacthes 
    mismatch_ch = get_chromosome.out.map{ dataset, vcf, chrms -> [ dataset, file(vcf), file(params.reference_genome) ] }
    check_mismatch(mismatch_ch)
    no_mismatchs = check_mismatch.out.map{ it -> no_mismatch(it) }

    // QC
    qc_ch = get_chromosome.out.map{ dataset, vcf, chrms -> [ dataset, file(chrms).readLines().unique().sort()[0], vcf ] }
    dataset_qc_dupl(qc_ch)
    split_multi_allelic(dataset_qc_dupl.out)
    filter_min_ac(split_multi_allelic.out.map{ dataset, chrm, vcf -> [ dataset, chrm, file(vcf), " --min-ac ${params.min_ac} --max-ac ${params.max_ac} "  ] })

    // Phasing without a reference
    phasing_ch = filter_min_ac.out.map{ dataset, chrm, vcf -> [ dataset, chrm, file(vcf), file(params.eagle_genetic_map) ] }
    phasing_vcf_no_ref(phasing_ch)
    vcf_to_m3vcf(phasing_vcf_no_ref.out)
    vcf_to_bcf(phasing_vcf_no_ref.out)
    vcf_legend(phasing_vcf_no_ref.out)



    // phasing_vcf(datasets_all_cha)
    // vcf_to_m3vcf(phasing_vcf.out)
    // vcf_to_bcf(phasing_vcf.out)
    // vcf_legend(phasing_vcf.out)

}