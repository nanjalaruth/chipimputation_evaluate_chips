#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process check_chromosome{

    input:
    val chrm

    output:
    tuple val(chrm), path("${chrm}.txt")

    """
    touch ${chrm} > ${chrm}.txt   
    """ 
}

workflow {
    chrms = Channel.from([1,22])
    check_chromosome(chrms)
    // Want to get a list to play with 
    check_chromosome.out.toList() // returns DataflowVariable(value=null) 
    check_chromosome.out.toList().getVal() // execution hangs up

}