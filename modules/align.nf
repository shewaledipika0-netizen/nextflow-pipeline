process ALIGN {

    publishDir "/home/payal/nextflow_pipeline/results/alignment", mode: 'copy'

    input:
    path reads
    path ref

    output:
    path "aligned.bam"

    script:
    """
    bwa index $ref
    bwa mem $ref $reads | samtools sort -o aligned.bam
    samtools index aligned.bam
    """
}
