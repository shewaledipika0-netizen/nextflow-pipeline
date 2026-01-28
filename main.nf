nextflow.enable.dsl=2

params.fastq = "/home/payal/nextflow_pipeline/data/fastq/*.fastq.gz"
params.ref   = "/home/payal/nextflow_pipeline/reference/GCF_000005845.2_ASM584v2_genomic.fna"
params.outdir = "/home/payal/nextflow_pipeline/results"

process FASTQC {
    tag "$reads"

    publishDir "results/fastqc", mode: 'symlink'


    input:
    path reads

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    fastqc ${reads}
    """
}



process ALIGNMENT {
    tag "$reads"

    input:
    path reads

    output:
    path "${reads.simpleName}.bam", emit: bam

    script:
    """
    bwa mem ${params.ref} ${reads} \
    | samtools view -Sb - \
    > ${reads.simpleName}.bam
    """
}

workflow {
    reads = Channel.fromPath(params.fastq)

    FASTQC(reads)
    ALIGNMENT(reads)
}
