process ALIGN {

    publishDir "/home/payal/nextflow_pipeline/results/alignment", mode: 'copy'

    input:
    path reads
    path ref

    output:
    path "aligned.bam"

    script:
    """
    ${params.bwa} index $ref
    ${params.bwa} mem $ref $reads | ${params.samtools} sort -o aligned.bam
    ${params.samtools} index aligned.bam
    """
}
