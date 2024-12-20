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
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
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
include { PHASE_GENOTYPES } from '../subworkflows/local/phase_genotypes'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: build locally
//
include { PYTHON3_FILTER_VCF } from '../modules/local/python3/filter_vcf/main'
include { PYTHON3_VCFTOLGEN } from '../modules/local/python3/vcftolgen/main'
include { PYTHON3_CREATE_GCTA_INPUT } from '../modules/local/python3/create_gcta_input/main'
include { CALC_UAR } from '../modules/local/uarprocess/calc_uar'
include { UAR_INVERSE } from '../modules/local/uarprocess/uar_inverse'
include { MAKE_GRM } from '../modules/local/gcta/make_grm'
include { REML_GRM } from '../modules/local/gcta/reml_grm'
include { H2_RANDOMPHENO_CREATE } from '../modules/local/randompheno/h2_randompheno_create'
include { REML_GRM as REML_GRM_SIM } from '../modules/local/gcta/reml_grm'
include { PLOT_H2_HISTOGRAM } from '../modules/local/gcta/plot_h2_histogram'
include { PYTHON3_PREPARE_PARAMETER_TEMPLATE } from '../modules/local/python3/prepare_parameter_template/main'
include { PYTHON3_CALC_LRT_ASREML } from '../modules/local/python3/calc_lrt_asreml/main'
include { PYTHON3_CALC_LRT_BLUPF90 } from '../modules/local/python3/calc_lrt_blupf90/main'
include { PYTHON3_PREPARE_MANHATTAN_PLOT_INPUT; PYTHON3_PREPARE_MANHATTAN_PLOT_INPUT as PROCESS_PERM } from '../modules/local/python3/prepare_manhttan_plot_input/main'
include { PYTHON3_MANHATTAN_PLOT } from '../modules/local/python3/manhattan_plot/main'
//
// MODULE: Installed directly from nf-core/modules
//
include { BCFTOOLS_CONCAT } from '../modules/nf-core/bcftools/concat/main'
include { PLINK2_VCF } from '../modules/nf-core/plink2/vcf/main'
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
    CHECK_INPUT(
        ch_input
    )

    if(!params.est_lrt_only){

    if (!params.include_chrom) {
        chrom = CHECK_INPUT.out.meta_vcf_idx.map { meta, vcf, idx -> meta.id } // chrom separated out so that UAR matrix can be calculated
    }

    else {
        include_chrom_file = Channel.fromPath(params.include_chrom)

        chrom = include_chrom_file.splitText().map { a->a.replaceAll("[\n\r]\$", '') }
    }

    pheno_f = Channel.fromPath(params.pheno_file) //read phenotype file

    if (params.skip_phasing == false) {
        PHASE_GENOTYPES(
        CHECK_INPUT.out.meta_vcf_idx
        )
        n1_meta_vcf_idx = PHASE_GENOTYPES.out.n1_meta_vcf_idx
    }
    else {
            n1_meta_vcf_idx = CHECK_INPUT.out.meta_vcf_idx
    }
    //
    // MODULE: filter vcf based on maf and sample inclusion-exclusion
    //

    PYTHON3_FILTER_VCF(
        n1_meta_vcf_idx.combine(pheno_f)
    )

    //PYTHON3_FILTER_VCF.out.gzvcf.view()

    vcf = PYTHON3_FILTER_VCF.out.gzvcf.map { meta, vcf -> vcf }

    outprefix = Channel.value(params.output_prefix)

    output_vcf = outprefix.combine(vcf)

    t_output_vcf = output_vcf.groupTuple() // tuple is grouped based on key of outprefix so that chromosome-splited vcf can be concatenated

    //
    // MODULE: concatenate chromosome-wise splited vcf file, key is output prefix
    //

    BCFTOOLS_CONCAT(
        t_output_vcf.map { prefix, vcf->tuple([id:prefix], vcf, []) }
    )

    // determine whether or not to estimate heritability, note that cldla and estimating heritability are mutually exclusive process

    if (params.estimate_h2) {
        PLINK2_VCF(
            BCFTOOLS_CONCAT.out.vcf
        )

        MAKE_GRM(
            PLINK2_VCF.out.bed,
            PLINK2_VCF.out.fam,
            PLINK2_VCF.out.bim
        )

        PYTHON3_CREATE_GCTA_INPUT(
            pheno_f
        )

        gi_gnb_gb = MAKE_GRM.out.gcta_grm_id.combine(MAKE_GRM.out.gcta_grm_n_bin).combine(MAKE_GRM.out.gcta_grm_bin)

        gi_gnb_gb_p = gi_gnb_gb.combine(PYTHON3_CREATE_GCTA_INPUT.out.gcta_phe)

        gi_gnb_gb_p_q_c = gi_gnb_gb_p.combine(PYTHON3_CREATE_GCTA_INPUT.out.gcta_qcovar.ifEmpty(null)).combine(PYTHON3_CREATE_GCTA_INPUT.out.gcta_covar.ifEmpty(null))

        REML_GRM(
            gi_gnb_gb_p_q_c.map { gi, gnb, gb, p, q, c->tuple(gi, gnb, gb, p, q == null ? [] : q, c == null ? [] : c) }
        )

        H2_RANDOMPHENO_CREATE(
            PYTHON3_CREATE_GCTA_INPUT.out.gcta_phe
        )

        ps = H2_RANDOMPHENO_CREATE.out.pheno_he.flatten()

        gi_gnb_gb_ps = gi_gnb_gb.combine(ps)

        gi_gnb_gb_ps_q_c = gi_gnb_gb_ps.combine(PYTHON3_CREATE_GCTA_INPUT.out.gcta_qcovar.ifEmpty(null)).combine(PYTHON3_CREATE_GCTA_INPUT.out.gcta_covar.ifEmpty(null))

        REML_GRM_SIM(
            gi_gnb_gb_ps_q_c.map { gi, gnb, gb, ps, q, c->tuple(gi, gnb, gb, ps, q == null ? [] : q, c == null ? [] : c) }
        )

        PLOT_H2_HISTOGRAM(
            REML_GRM.out.gcta_hsq,
            REML_GRM_SIM.out.gcta_hsq.collect()
        )
    }

    else {
        mergedvcf_chrom = BCFTOOLS_CONCAT.out.vcf.map { chrom, vcf->vcf }.combine(chrom) // use "chrom" value splitted on line number 81

        chrom_mergedvcf = mergedvcf_chrom.map { vcf, chrom -> tuple(chrom, vcf) } // rearrange tuple with chrom as key

        //
        // MODULE: convert vcf to lgen so that lgen file can be supplied to calculate uar matrix
        //
        PYTHON3_VCFTOLGEN(
        chrom_mergedvcf.combine(pheno_f)
    )

    //
    // MODULE: calculate UAR matrix chromosome-wise for each lgen file generated in previous module
    //

        CALC_UAR(
        PYTHON3_VCFTOLGEN.out.chrom_lgen
    )

        chrom_lgen_uar = PYTHON3_VCFTOLGEN.out.chrom_lgen.combine(CALC_UAR.out.chrom_uar, by:0) //lgen file is supplied --> inverse require total number of indi

        //
        // MODULE: calculate inverse of UAR matrix
        //
        UAR_INVERSE(
        chrom_lgen_uar
    )

    //
    // MODULE: prepare parameter template file for the tools asreml or blupf90
    //

        par_file = params.extra_params ? Channel.fromPath(params.extra_params) : (params.tool == 'asreml' ? Channel.fromPath(params.default_asreml_params) : Channel.fromPath(params.default_blupf90_params))

        PYTHON3_PREPARE_PARAMETER_TEMPLATE(
        pheno_f,
        par_file
    )

        ch_lrt = UAR_INVERSE.out.giv.combine(PYTHON3_PREPARE_PARAMETER_TEMPLATE.out.updated_phe).combine(params.tool == 'asreml' ? PYTHON3_PREPARE_PARAMETER_TEMPLATE.out.asr_template : PYTHON3_PREPARE_PARAMETER_TEMPLATE.out.blp_template).combine(PYTHON3_FILTER_VCF.out.gzvcf,by : 0)

    //
    // MODULE: main python script to calculate lrt values
    //

        PYTHON3_CALC_LRT = params.tool == 'asreml' ? PYTHON3_CALC_LRT_ASREML(ch_lrt) : PYTHON3_CALC_LRT_BLUPF90(ch_lrt)


        lrt_f = PYTHON3_CALC_LRT.real_txt.map { v, f->f }.collect()
        lrt_perm_f = PYTHON3_CALC_LRT.perm_txt.map { v, f->f }.collect()

        //
        // MODULE: python script to prepare the input file for manhattan plot
        //
        PYTHON3_PREPARE_MANHATTAN_PLOT_INPUT(
        ch_input,
        lrt_f
    )

        //
        // MODULE: python script to prepare the input file for manhattan plot
        //
        PROCESS_PERM(
        ch_input,
        lrt_perm_f
    )

        plot_yml = Channel.fromPath(params.manhattan_plot_yml, checkIfExists: true)

    //
    //MODULE: for generating manhattan plot
    //

        plot_yml = Channel.fromPath(params.manhattan_plot_yml, checkIfExists: true)

        PYTHON3_MANHATTAN_PLOT(
        PYTHON3_PREPARE_MANHATTAN_PLOT_INPUT.out.maninp_txt,
        plot_yml,
        PROCESS_PERM.out.threshold_txt
    )
        }
    }
    // this "else" loop will only run when "--est-lrt-only" is set to true, decided not to touch the modules and change the channel output of the part of the worklfow that are already running smoothly
    else{

    //
    // MODULE: prepare parameter template file for the tools asreml or blupf90
    //

        par_file = params.extra_params ? Channel.fromPath(params.extra_params) : (params.tool == 'asreml' ? Channel.fromPath(params.default_asreml_params) : Channel.fromPath(params.default_blupf90_params))

        PYTHON3_PREPARE_PARAMETER_TEMPLATE(
        pheno_f,
        par_file
        )

        }
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
