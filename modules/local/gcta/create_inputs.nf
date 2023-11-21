process CREATE_INPUTS{

    tag { "creating inputs for GCTA" }
    label "process_medium"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/gcta/inputs/", mode:"copy")

    input:
        path(pheno_f)
        path(qcovar_f)
        path(covar_f)

    output:
        path ("*.phe"), emit: gcta_phe
        path ("*.qcovar"), emit: gcta_qcovar, optional: true
        path ("*.covar"), emit: gcta_covar, optional: true
        
    
    script:
        def pheno_col = params.col_pheno
        output_prefix = params.output_prefix

        if (qcovar_f != [] && covar_f != []){

        """
        awk -v p_col=${pheno_col} '{print \$2,\$2,\$p_col}' ${pheno_f} > ${output_prefix}_gcta_in.phe
        
        awk 'NR==FNR{a[\$0];next}{printf \$2" "\$2;for(i in a){printf " "\$i};printf "\\n"}' ${qcovar_f} ${pheno_f} > ${output_prefix}_gcta_in.qcovar

        awk 'NR==FNR{a[\$0];next}{printf \$2" "\$2;for(i in a){printf " "\$i};printf "\\n"}' ${covar_f} ${pheno_f} > ${output_prefix}_gcta_in.covar
        

        """ 
        }
        else{
            if( qcovar_f !=[] && covar_f == []){

            """

            awk -v p_col=${pheno_col} '{print \$2,\$2,\$p_col}' ${pheno_f} > ${output_prefix}_gcta_in.phe
            
            awk 'NR==FNR{a[\$0];next}{printf \$2" "\$2;for(i in a){printf " "\$i};printf "\\n"}' ${qcovar_f} ${pheno_f} > ${output_prefix}_gcta_in.qcovar

            """


            }

            else{
                if ( qcovar_f == [] && covar_f != [] ){

                """

        awk -v p_col=${pheno_col} '{print \$2,\$2,\$p_col}' ${pheno_f} > ${output_prefix}_gcta_in.phe
        
        awk 'NR==FNR{a[\$0];next}{printf \$2" "\$2;for(i in a){printf " "\$i};printf "\\n"}' ${covar_f} ${pheno_f} > ${output_prefix}_gcta_in.covar


                """


                }

                else{
                """

                    awk -v p_col=${pheno_col} '{print \$2,\$2,\$p_col}' ${pheno_f} > ${output_prefix}_gcta_in.phe

                """


                }


            }

        }
}
