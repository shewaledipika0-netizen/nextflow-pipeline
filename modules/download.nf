process DOWNLOAD {

    publishDir "/home/payal/nextflow_pipeline/data", mode: 'copy'

    output:
    path "reads.fastq"
    path "reference.fa"

    script:
    """
    # create valid FASTQ file
    echo "@read1
ACTGACTGACTGACTGACTGACTG
+
FFFFFFFFFFFFFFFFFFFFFFFF" > reads.fastq

    # matching reference genome
    echo ">chr1
ACTGACTGACTGACTGACTGACTG" > reference.fa
    """
}
