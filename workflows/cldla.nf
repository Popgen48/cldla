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
include { CALC_LRT } from '../subworkflows/local/calc_lrt'
include { CALC_LRT as CALC_LRT_PERMUTATION } from '../subworkflows/local/calc_lrt'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: build locally 
//

include { VCFTOLGEN } from '../modules/local/vcftolgen/vcftolgen'
include { CALC_UAR } from '../modules/local/uarprocess/calc_uar'
include { UAR_INVERSE } from '../modules/local/uarprocess/uar_inverse'
include { CREATEHAPLO_RECORD } from '../modules/local/createhaplo/createhaplo_record'
include { CREATE_HAPMAP } from '../modules/local/createhaplo/create_hapmap'
include { CREATE_INVERSE } from '../modules/local/grmprocess/create_inverse'
include { RANDOMPHENO_CREATE } from '../modules/local/randompheno/randompheno_create'
include { MAKE_GRM } from '../modules/local/gcta/make_grm'
include { REML_GRM } from '../modules/local/gcta/reml_grm'
include { CREATE_INPUTS } from '../modules/local/gcta/create_inputs'
include { H2_RANDOMPHENO_CREATE } from '../modules/local/randompheno/h2_randompheno_create'
include { REML_GRM as REML_GRM_SIM } from '../modules/local/gcta/reml_grm'
include { PLOT_H2_HISTOGRAM } from '../modules/local/gcta/plot_h2_histogram'
//
// MODULE: Installed directly from nf-core/modules
//
include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main'
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

        chrom = CHECK_INPUT.out.chrom_vcf_idx.map{chrom, vcf, idx -> chrom} // chrom separated out so that UAR matrix can be calculated 

        }

    else{

        include_chrom_file = Channel.fromPath(params.include_chrom)

        chrom = include_chrom_file.splitText().map{a->a.replaceAll("[\n\r]\$", "")}

        }
        
        
    pheno_f = Channel.fromPath(params.pheno_file) //read phenotype file

    //
    // MODULE: filter vcf based on maf and sample inclusion-exclusion
    //
    
    //CHECK_INPUT.out.chrom_vcf_idx.view()

    BCFTOOLS_VIEW(
        CHECK_INPUT.out.chrom_vcf_idx,
        [],
        [],
        []
    )


    vcf = BCFTOOLS_VIEW.out.vcf.map{chrom,vcf -> vcf}

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
    VCFTOLGEN(
        chrom_mergedvcf
    )
    
    //
    // MODULE: calculate UAR matrix chromosome-wise for each lgen file generated in previous module
    //

    CALC_UAR(
        VCFTOLGEN.out.chrom_lgen
    )
    
    chrom_lgen_uar = VCFTOLGEN.out.chrom_lgen.combine(CALC_UAR.out.chrom_uar, by:0) //lgen file is supplied --> inverse require total number of indi

    //
    // MODULE: calculate inverse of UAR matrix
    //
    UAR_INVERSE(
        chrom_lgen_uar
    )

    //
    // MODULE: create total windows to be generated for each chromosome given window size
    //
    CREATEHAPLO_RECORD(
        BCFTOOLS_VIEW.out.vcf.combine(chrom,by:0)
    )
    chrom_window_out = CREATEHAPLO_RECORD.out.record.splitText().map{a->tuple(a.split()[0], a.split()[1], a.split()[2])} //read window file line by line

    chrom_window_out_vcf = chrom_window_out.combine(BCFTOOLS_VIEW.out.vcf, by: 0) // combine filtered vcf file 


    //
    // MODULE: create hap, map and par file for each window each chromosome given vcf
    //

    CREATE_HAPMAP(
        chrom_window_out_vcf
    )


    //
    // MODULE: calculate grm, its bending and get inverse of grm for each window given outputs from the previous module
    //

    CREATE_INVERSE(
        CREATE_HAPMAP.out.chrom_hap,
        CREATE_HAPMAP.out.chrom_map,
        CREATE_HAPMAP.out.chrom_par
    )

    chrwin_hap = CREATE_HAPMAP.out.chrom_hap.map{chrom,hap->tuple(hap.getBaseName().minus('.Hap'), hap)} //rearrange key of hap channel
    
    chrwin_map = CREATE_HAPMAP.out.chrom_map.map{chrom,map_f->tuple(map_f.getBaseName().minus('.Map'), map_f)} //rearrange key of map channel

    chrwin_hap_map = chrwin_hap.combine(chrwin_map, by:0) // combine hap and map channel of the same chromosome, same window

    chrwin_hap_map_winginv = chrwin_hap_map.combine(CREATE_INVERSE.out.chrwin_ginv, by:0) //combine ginv, hap and map of the same chromosome same window

    chr_hap_map_winginv = chrwin_hap_map_winginv.map{chrwin, hap, map_f, ginv->tuple(chrwin.split("\\.")[-2], hap, map_f, ginv)} // rearrange key to chromosome id

    chr_hap_map_winginv_chromginv = chr_hap_map_winginv.combine(UAR_INVERSE.out.chrom_ginv,by:0) // combine inverse of uar matrix of the same chromosome (key is chromosome)

    chr_hap_map_winginv_chromginv_pheno = chr_hap_map_winginv_chromginv.combine(pheno_f)

    cldla_val = Channel.value("cldla")

    CALC_LRT(
        chr_hap_map_winginv_chromginv_pheno,
        cldla_val
        )

    }
    if ( params.permutation_test && !params.estimate_h2 ){
                //read file of random chromosome window generated by CREATEHAPLO_RECORD module, important to remove newline character to match key
                chrmwin_outprefix = CREATEHAPLO_RECORD.out.random_record.splitText().map{a->a.replaceAll("[\n\r]\$", "")}.combine(outprefix)

                out_chromwin = chrmwin_outprefix.map{chromwin, outprefix->tuple(outprefix+"."+chromwin, chromwin)} // rearrange key

                out_chromwin_ginv = out_chromwin.combine(CREATE_INVERSE.out.chrwin_ginv, by : 0) // combine ginv of matching key of random window

                out_chromwin_ginv_hap_map = out_chromwin_ginv.combine(chrwin_hap_map,by:0) // combine .hap and .map file of the same random window

                chrom_ginv_hap_map = out_chromwin_ginv_hap_map.map{out, chromwin, ginv, hap, map_f->tuple(chromwin.split("\\.")[-2],ginv,hap,map_f)} //make chrm as key

                chrom_pheno = chrom.combine(pheno_f) //use chrom value of line 81

                //
                // MODULE: to generate "n" datasets of randomized phenotypes for the permutation test
                //

                RANDOMPHENO_CREATE(
                    chrom_pheno
                )

                //flatten all elements of tuple, and make chromosome as key 
                chrom_ranpheno = RANDOMPHENO_CREATE.out.pheno_r.flatten().map{i->tuple(i.getBaseName().split("_")[-2], i)} 

                chrom_ginv_hap_map_ranpheno = chrom_ginv_hap_map.combine(chrom_ranpheno, by:0) //merge data channel with channel of chrom, winginv, hap

                chrom_ginv_hap_map_ranpheno_uarinv = chrom_ginv_hap_map_ranpheno.combine(UAR_INVERSE.out.chrom_ginv,by:0)

                permut_val = Channel.value('permute')

                CALC_LRT_PERMUTATION(
                    chrom_ginv_hap_map_ranpheno_uarinv.map{chrom, ginv, hap, map_f, phe, uar->tuple(chrom, hap, map_f, ginv, uar, phe)},
                    permut_val
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
