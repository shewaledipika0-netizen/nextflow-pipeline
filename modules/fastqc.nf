process FASTQC {

    publishDir "/home/payal/nextflow_pipeline/results/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "*.html"

    script:
    """
    fastqc $reads
    """
}
