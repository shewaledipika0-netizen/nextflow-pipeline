process VARIANT {

    publishDir "/home/payal/nextflow_pipeline/results/variants", mode: 'copy'

    input:
    path bam
    path ref

    output:
    path "variants.vcf"

    script:
    """
    bcftools mpileup -f $ref $bam | bcftools call -mv -Ov -o variants.vcf
    """
}
