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
// MODULE: build locally 
//
include { PYTHON3_FILTER_VCF } from '../modules/local/python3/filter_vcf/main'
include { PYTHON3_VCFTOLGEN } from '../modules/local/python3/vcftolgen/main'
include { CALC_UAR } from '../modules/local/uarprocess/calc_uar'
include { UAR_INVERSE } from '../modules/local/uarprocess/uar_inverse'
include { MAKE_GRM } from '../modules/local/gcta/make_grm'
include { REML_GRM } from '../modules/local/gcta/reml_grm'
include { CREATE_INPUTS } from '../modules/local/gcta/create_inputs'
include { H2_RANDOMPHENO_CREATE } from '../modules/local/randompheno/h2_randompheno_create'
include { REML_GRM as REML_GRM_SIM } from '../modules/local/gcta/reml_grm'
include { PLOT_H2_HISTOGRAM } from '../modules/local/gcta/plot_h2_histogram'
include { PYTHON3_CALC_LRT } from '../modules/local/python3/calc_lrt/main'
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
    CHECK_INPUT (
        ch_input
    )

    if( params.include_chrom == "none"){

        chrom = CHECK_INPUT.out.meta_vcf_idx.map{meta, vcf, idx -> meta.id} // chrom separated out so that UAR matrix can be calculated 

        }

    else{

        include_chrom_file = Channel.fromPath(params.include_chrom)

        chrom = include_chrom_file.splitText().map{a->a.replaceAll("[\n\r]\$", "")}

        }
        
        
    pheno_f = Channel.fromPath(params.pheno_file) //read phenotype file


    //
    // MODULE: filter vcf based on maf and sample inclusion-exclusion
    //

    PYTHON3_FILTER_VCF(
        CHECK_INPUT.out.meta_vcf_idx.combine(pheno_f)
    )

    //PYTHON3_FILTER_VCF.out.gzvcf.view()
    


    vcf = PYTHON3_FILTER_VCF.out.gzvcf.map{meta,vcf -> vcf}

    outprefix = Channel.value(params.output_prefix)

    output_vcf = outprefix.combine(vcf)

    t_output_vcf = output_vcf.groupTuple() // tuple is grouped based on key of outprefix so that chromosome-splited vcf can be concatenated

    //
    // MODULE: concatenate chromosome-wise splited vcf file, key is output prefix
    //

    BCFTOOLS_CONCAT(
        t_output_vcf.map{prefix,vcf->tuple(prefix,vcf,[])}
    )

    // determine whether or not to estimate heritability, note that cldla and estimating heritability are mutually exclusive process

    if ( params.estimate_h2 ){
        

        PLINK2_VCF(
            BCFTOOLS_CONCAT.out.vcf
        )

        MAKE_GRM(
            PLINK2_VCF.out.bed,
            PLINK2_VCF.out.fam,
            PLINK2_VCF.out.bim
        )
        
    

        qcovar_file = Channel.fromPath(params.qcovar)
        
        covar_file = Channel.fromPath(params.covar)

        CREATE_INPUTS(
            pheno_f,
            qcovar_file.map{q->q=="none"?[]:q},
            covar_file.map{c->c=="none"?[]:c}
        )

        gi_gnb_gb = MAKE_GRM.out.gcta_grm_id.combine(MAKE_GRM.out.gcta_grm_n_bin).combine(MAKE_GRM.out.gcta_grm_bin)

        gi_gnb_gb_p_q_c = gi_gnb_gb.combine(CREATE_INPUTS.out.gcta_phe).combine(CREATE_INPUTS.out.gcta_qcovar).combine(CREATE_INPUTS.out.gcta_covar)
        

        REML_GRM(
            gi_gnb_gb_p_q_c.map{gi,gnb,gb,p,q,c->tuple(gi,gnb,gb,p,q=="none"?[]:q,c=="none"?[]:c)}
        )

        H2_RANDOMPHENO_CREATE(
            CREATE_INPUTS.out.gcta_phe
        )
        
        ps = H2_RANDOMPHENO_CREATE.out.pheno_he.flatten()

        gi_gnb_gb_ps_q_c = gi_gnb_gb.combine(ps).combine(CREATE_INPUTS.out.gcta_qcovar).combine(CREATE_INPUTS.out.gcta_covar)

        REML_GRM_SIM(
            gi_gnb_gb_ps_q_c.map{gi,gnb,gb,ps,q,c->tuple(gi,gnb,gb,ps,q=="none"?[]:q,c=="none"?[]:c)}
        )

        PLOT_H2_HISTOGRAM(
            REML_GRM.out.gcta_hsq,
            REML_GRM_SIM.out.gcta_hsq.collect()
        )
        

    }

    else{

    mergedvcf_chrom = BCFTOOLS_CONCAT.out.vcf.map{chrom,vcf->vcf}.combine(chrom) // use "chrom" value splitted on line number 81

    chrom_mergedvcf = mergedvcf_chrom.map{vcf, chrom -> tuple(chrom, vcf)} // rearrange tuple with chrom as key

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

    par_file = Channel.fromPath(params.par_file) //read parameter file

    ch_lrt = UAR_INVERSE.out.giv.combine(pheno_f).combine(par_file).combine(PYTHON3_FILTER_VCF.out.gzvcf,by:0)

    PYTHON3_CALC_LRT(
        ch_lrt
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
