include { RUN_ECHIDNA as RUN_ECHIDNA_H1 } from '../../modules/local/echidna/run_echidna'
include { RUN_ECHIDNA as RUN_ECHIDNA_H0 } from '../../modules/local/echidna/run_echidna'
include { RUN_ECHIDNA_SERIAL as RUN_ECHIDNA_SERIAL_H1 } from '../../modules/local/echidna/run_echidna_serial'
include { RUN_ECHIDNA as RUN_ECHIDNA_PERMUTATION } from '../../modules/local/echidna/run_echidna'
include { RUN_ASREML_SERIAL } from '../../modules/local/asreml/run_asreml_serial'
include { RUN_ASREML_WINDOW as RUN_ASREML_H1 } from '../../modules/local/asreml/run_asreml_window.nf'
include { RUN_ASREML_WINDOW as RUN_ASREML_H0 } from '../../modules/local/asreml/run_asreml_window.nf'
include { RUN_ASREML_SERIAL as RUN_ASREML_SERIAL_PERMUTE } from '../../modules/local/asreml/run_asreml_serial'
include { CALC_LRT_RATIO } from '../../modules/local/calc_lrt_ratio'
include { PLOT_LRT_RATIO } from '../../modules/local/plot_lrt_ratio'
include { RUN_AIREMLF90_SERIAL } from '../../modules/local/airemlf90/run_airemlf90_serial'
include { RUN_AIREMLF90_WINDOW_WISE as RUN_AIREMLF90_H0 } from '../../modules/local/airemlf90/run_airemlf90_window_wise'
include { RUN_AIREMLF90_WINDOW_WISE as RUN_AIREMLF90_H1 } from '../../modules/local/airemlf90/run_airemlf90_window_wise'

workflow CALC_LRT{
    take:
        chr_hap_map_winginv_chromginv_pheno
        w_value

    main:
        if( params.run_window_serial ){

            chr_chrinv_phe = chr_hap_map_winginv_chromginv_pheno.map{chr,hap,map,winginv,chromginv,pheno->tuple(chr,chromginv,pheno)}.unique()

            if ( params.echidna ){
            //
            // MODULE: run echidna for H1 hypothesis i.e. for each window
            //
            RUN_ECHIDNA_H1(
                chr_hap_map_winginv_chromginv_pheno
            )

            //
            // MODULE: run echidna for H0 hypothesis i.e. for each chromosome without ginv file of window 
            //
        
            RUN_ECHIDNA_H0(
                chr_chrinv_phe.map{chr, chrinv, phe->tuple(chr, [], [], [], chrinv, phe)}
            )

            chrom_h0llk_h1llk = RUN_ECHIDNA_H0.out.chrom_llik.combine(RUN_ECHIDNA_H1.out.chrom_llik, by:0)
            
            CALC_LRT_RATIO(
                chrom_h0llk_h1llk
            )
            }
        
            if ( params.airemlf90 ){

            //

            // MODULE: run echidna for H1 hypothesis i.e. for each window
            //
            RUN_AIREMLF90_H1(
                chr_hap_map_winginv_chromginv_pheno,
                w_value
            )

            //
            // MODULE: run echidna for H0 hypothesis i.e. for each chromosome without ginv file of window 
            //
        
            RUN_AIREMLF90_H0(
                chr_chrinv_phe.map{chr, chrinv, phe->tuple(chr, [], [], [], chrinv, phe)},
                w_value
            )

            CALC_LRT_RATIO(
                RUN_AIREMLF90_H0.out.chrom_llik.combine(RUN_AIREMLF90_H1.out.chrom_llik.groupTuple(), by:0)
            )

            PLOT_LRT_RATIO(
                CALC_LRT_RATIO.out.lrt.collect()
            )


            }

            if( params.asreml ) {

            //
            // MODULE: run echidna for H1 hypothesis i.e. for each window
            //
            RUN_ASREML_H1(
                chr_hap_map_winginv_chromginv_pheno
            )

            //
            // MODULE: run echidna for H0 hypothesis i.e. for each chromosome without ginv file of window 
            //
        
            RUN_ASREML_H0(
                chr_chrinv_phe.map{chr, chrinv, phe->tuple(chr, [], [], [], chrinv, phe)}
            )

            chrom_h0llk_h1llk = RUN_ASREML_H0.out.chrom_llik.combine(RUN_ASREML_H1.out.chrom_llik.groupTuple(), by:0)
            
            CALC_LRT_RATIO(
                chrom_h0llk_h1llk
            )


            }

        }
        else{
            
            c_h_m_w_c_p = chr_hap_map_winginv_chromginv_pheno.groupTuple().map{c,h,m,w,ci,pheno->tuple(c,h.unique(),m.unique(),w.unique(),ci.unique(),pheno)}

        if( params.echidna ){
    
            RUN_ECHIDNA_SERIAL_H1(
                c_h_m_w_c_p,
                w_value
            )
        }

        if( params.asreml ){
                    //c_h_m_w_c_pi = chr_hap_map_winginv_chromginv_pheno.groupTuple().map{c,h,m,w,ci,pheno->tuple(c,h.unique(),m.unique(),w.unique(),ci.unique(),pheno.unique())}
                    RUN_ASREML_SERIAL(
                        c_h_m_w_c_p,
                        w_value
                        )
                    PLOT_LRT_RATIO(
                        RUN_ASREML_SERIAL.out.lrt_out.collect()
                    )
                        
                }
        if( params.airemlf90 ){

                    //c_h_m_w_c_pi = chr_hap_map_winginv_chromginv_pheno.groupTuple().map{c,h,m,w,ci,pheno->tuple(c,h.unique(),m.unique(),w.unique(),ci.unique(),pheno.unique())}
                    RUN_AIREMLF90_SERIAL(
                        c_h_m_w_c_p,
                        w_value
                        )
        }
    }
                

}
