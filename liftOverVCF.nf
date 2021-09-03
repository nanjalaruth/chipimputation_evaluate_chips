#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}


process sort_vcf{
    tag "sort_vcf_${dataset}_${chrm}"
    label "bcftools"

    input:
      tuple val(dataset), file(vcf_file), val(chrm)

    output:
      tuple val(dataset), file(vcf_out), val(chrm)

    script:
      vcf_out = "${vcf_file.getSimpleName()}.sorted.bcf"
      """
      bcftools view ${vcf_file} -Ob -o ${vcf_file.getSimpleName()}.bcf
      bcftools sort ${vcf_file.getSimpleName()}.bcf -T . -Ob -o ${vcf_out}
      """
}

process crossmap{
    tag "crossmap_${dataset}_${chrm}"
    label "crossmap"

    input:
      tuple val(dataset), file(vcf_file), val(chrm), file(chain_file), file(reference_genome)

    output:
      tuple val(dataset), file(vcf_out), val(chrm)

    script:
      vcf_out = "${vcf_file.getSimpleName()}.crossmap.vcf"
      """
      CrossMap.py vcf ${chain_file} ${vcf_file} ${reference_genome} ${vcf_out}
      """
}

process convert_to_zero_based_bed{
    tag "to_zero_bed_${dataset}"
    label "crossmap"

    input:
      tuple val(dataset), file(bed_file)

    output:
      tuple val(dataset), file(bed_out)

    script:
      bed_out = "${bed_file.getSimpleName()}.zero-based.b38.bed"
      """
      awk '{print \$1"\\t"\$2-1"\\t"\$3}' ${bed_file} > ${bed_out}
      """
}

process crossmap_bed_files{
    tag "crossmap_bed_${dataset}"
    label "crossmap"

    input:
      tuple val(dataset), file(bed_file), file(chain_file)

    output:
      tuple val(dataset), file(bed_out)

    script:
      bed_out = "${bed_file.getSimpleName()}.crossmap.b38.bed"
      """
      CrossMap.py bed ${chain_file} ${bed_file} ${bed_out}
      """
}

process rename_chrs{
    tag "rename-chrs_${dataset}_${chrm}"
    label "bcftools"

    input:
      tuple val(dataset), file(vcf_file), val(chrm), file(chromosome_mapping_file)

    output:
      tuple val(dataset), file(vcf_out), val(chrm)

    script:
      vcf_out = "${vcf_file.getSimpleName()}.rename-chrs.bcf"
      """
      bcftools annotate --rename-chrs ${chromosome_mapping_file} ${vcf_file} -Ob -o ${vcf_out}
      """
}

process add_chr_bed{
    tag "add_chr_bed_${dataset}"
    label "crossmap"
    publishDir "${params.outdir}/chip_files", overwrite: true, mode:'copy', pattern: "*.bed"

    input:
      tuple val(dataset), file(bed_file)

    output:
      tuple val(dataset), file(bed_out)

    script:
      bed_out = "${bed_file.getSimpleName()}.b38.bed"
      """
      awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' ${bed_file} > ${bed_out}
      """
}


process fix_contig{
    tag "fix_contig_${dataset}_${chrm}"
    label "bcftools"
    publishDir "${params.outdir}/b38", overwrite: false, mode:'copy', pattern: "*vcf.gz*"

    input:
      tuple val(dataset), file(vcf_file), val(chrm), file(reference_genome)

    output:
      tuple val(dataset), file(vcf_out), val(chrm)

    script:
      vcf_out = "${vcf_file.getSimpleName()}.b38.vcf.gz"
      """
      samtools faidx ${reference_genome}
      bcftools view ${vcf_file} -o ${vcf_file.getSimpleName()}.vcf
      bcftools reheader -f ${reference_genome}.fai ${vcf_file.getSimpleName()}.vcf -o ${file(vcf_out).getSimpleName()}.temp.vcf
      bcftools view ${file(vcf_out).getSimpleName()}.temp.vcf -Oz -o ${vcf_out}
      rm ${file(vcf_out).getSimpleName()}*vcf
      """
}

process vcf_index{
    tag "index_${dataset}_${chrm}"
    label "bcftools"
    publishDir "${params.outdir}/b38", overwrite: false, mode:'copy', pattern: "*vcf.gz.tbi*"

    input:
      tuple val(dataset), file(vcf_file), val(chrm)

    output:
      tuple val(dataset), file(vcf_file), file("${vcf_file}.tbi"), val(chrm)

    script:
      """
      tabix -f ${vcf_file}
      """
}

// workflow vcf_files{

// }

workflow bed_files{
  take:
    data
  main:
    convert_to_zero_based_bed( data )
    crossmap_bed_files( convert_to_zero_based_bed.out.map{ dataset, bed_file -> [ dataset, file(bed_file), file(params.chain_file) ]} )
    add_chr_bed( crossmap_bed_files.out )
  emit:
    data

}

workflow{

  //// For VCF files
  datasets_all = []
  params.vcf_files.each { dataset, dataset_vcf ->
      check_files([dataset_vcf])
      datasets_all << [ dataset, file(dataset_vcf) ]
  }
  datasets_all_cha = Channel.from(datasets_all)
  crossmap( datasets_all_cha.map{ dataset, vcf_file -> [ dataset, file(vcf_file), '', file(params.chain_file), file(params.reference_genome) ]} )
  sort_vcf( crossmap.out.map{ dataset, vcf_file, chrm -> [ dataset, file(vcf_file), '' ] } )
  rename_chrs( sort_vcf.out.map{ dataset, vcf_file, chrm -> [ dataset, file(vcf_file), '', file(params.chromosome_mapping)]  } )
  fix_contig( rename_chrs.out.map{ dataset, vcf_file, chrm -> [ dataset, file(vcf_file), '', file(params.reference_genome) ] } )
  vcf_index( fix_contig.out )

  //// For MAP files or o-based BED files ----> chr start end
  map_files = []
  params.map_files.each { dataset, dataset_bed ->
      check_files([dataset_bed, params.chain_file])
      map_files << [ dataset, file(dataset_bed) ]
  }
  map_files_cha = Channel.from(map_files)
  bed_files( map_files_cha )

}

