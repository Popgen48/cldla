/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCldla.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { CHECK_INPUT } from '../subworkflows/local/check_input'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_CONCAT } from '../modules/nf-core/bcftools/concat/main'
include { VCFTOLGEN } from '../modules/local/vcftolgen/vcftolgen'
include { CALC_UAR } from '../modules/local/uarprocess/calc_uar'
include { UAR_INVERSE } from '../modules/local/uarprocess/uar_inverse'
include { CREATEHAPLO_RECORD } from '../modules/local/createhaplo/createhaplo_record'
include { CREATE_HAPMAP } from '../modules/local/createhaplo/create_hapmap'
include { CREATE_INVERSE } from '../modules/local/grmprocess/create_inverse'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CLDLA {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    CHECK_INPUT (
        ch_input
    )

    chrom = CHECK_INPUT.out.chrom_vcf_idx.map{chrom, vcf, idx -> chrom}

    BCFTOOLS_VIEW(
        CHECK_INPUT.out.chrom_vcf_idx,
        [],
        [],
        []
    )

    vcf = BCFTOOLS_VIEW.out.vcf.map{chrom,vcf -> vcf}
    outprefix = Channel.value(params.output_prefix)
    output_vcf = outprefix.combine(vcf)
    t_output_vcf = output_vcf.groupTuple()
    t_output_vcf.view()

    BCFTOOLS_CONCAT(
        t_output_vcf.map{prefix,vcf->tuple(prefix,vcf,[])}
    )

    mergedvcf_chrom = BCFTOOLS_CONCAT.out.vcf.map{chrom,vcf->vcf}.combine(chrom)

    chrom_mergedvcf = mergedvcf_chrom.map{vcf, chrom -> tuple(chrom, vcf)}

    //chrom_mergedvcf.view()

    VCFTOLGEN(
        chrom_mergedvcf
    )

    CALC_UAR(
        VCFTOLGEN.out.chrom_lgen
    )
    
    chrom_lgen_uar = VCFTOLGEN.out.chrom_lgen.combine(CALC_UAR.out.chrom_uar, by:0)

    //chrom_lgen_uar.view()

    UAR_INVERSE(
        chrom_lgen_uar
    )

    //BCFTOOLS_VIEW.out.vcf.view()

    CREATEHAPLO_RECORD(
        BCFTOOLS_VIEW.out.vcf
    )
    chrom_window = CREATEHAPLO_RECORD.out.record.splitText().map{a->tuple(a.split()[0], a.split()[1])}

    chrom_window_vcf = chrom_window.combine(BCFTOOLS_VIEW.out.vcf, by: 0)

    CREATE_HAPMAP(
        chrom_window_vcf
    )
    //CREATE_HAPMAP.out.record.view()
    //chrom_hap_map_par = CREATE_HAPMAP.out.record.map{chrom,hmp->tuple(chrom,hmp[0],hmp[1],hmp[2])}
    //chrom_hap_map_par.view()
    CREATE_INVERSE(
        CREATE_HAPMAP.out.chrom_hap,
        CREATE_HAPMAP.out.chrom_map,
        CREATE_HAPMAP.out.chrom_par
    )

    chrwin_hap = CREATE_HAPMAP.out.chrom_hap.map{chrom,hap->tuple(hap.getBaseName().minus('.Hap'), hap)}
    chrwin_hap_ginv = chrwin_hap.combine(CREATE_INVERSE.out.chrwin_ginv, by:0)
    chrwin_hap_ginv.view()
    //chrom_windowinv_chrominv = UAR_INVERSE.out.chrom_ginv.combine(CREATE_INVERSE.out.chrom_ginv, by:0)
    //chrom_windowinv_chrominv_hap = chrom_windowinv_chrominv.combine(CREATE_HAPMAP.out.chrom_hap, by:0)
    //chrom_windowinv_chrominv_hap.view()
    
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
