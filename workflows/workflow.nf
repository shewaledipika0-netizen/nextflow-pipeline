include { DOWNLOAD } from '../modules/download.nf'
include { FASTQC } from '../modules/fastqc.nf'
include { ALIGN } from '../modules/align.nf'
include { VARIANT } from '../modules/variant.nf'

workflow VARIANT_WORKFLOW {

    files = DOWNLOAD()

    FASTQC(files[0])
    aligned = ALIGN(files[0], files[1])
    VARIANT(aligned, files[1])
}
