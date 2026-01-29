nextflow.enable.dsl=2

include { VARIANT_WORKFLOW } from './workflows/workflow.nf'

workflow {
    VARIANT_WORKFLOW()
}
