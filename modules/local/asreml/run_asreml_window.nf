process RUN_ASREML_WINDOW{

    tag { "running_echidna_${chrom}" }
    label "process_medium"
    maxForks 1
    publishDir("${params.outdir}/echidna/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(hap), path(map_f), path(win_ginv), path(chrom_ginv), path(pheno_file)

    output:
        tuple val("${chrom}"), path ("*.{as,asr}"), emit: chrom_as
        tuple val("${chrom}"), path ( "*.llik" ), emit: chrom_llik
        
    
    script:
        
        echidna_params = params.echidna_params
        win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        chrom_ginv_base = chrom_ginv.getName()
        new_prefix = hap==[] ? chrom_ginv.getName().minus('.giv'):hap.getName().minus('.Hap')
        window_size = params.window_size

        if( hap != [] ){
        
        """

        diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' ${hap})

        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' ${hap} ${pheno_file} > ${new_prefix}.phe
        
        python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} \${diplo} ${chrom_ginv_base} ${win_ginv_base} ${new_prefix}.phe h1

        asreml -NS5 ${new_prefix}.as

        python ${baseDir}/bin/extract_loglik.py ${new_prefix}.asr ${map_f} ${new_prefix}.llik ${chrom} ${window_size} asreml


        """ 
        }
        else{

        """

        python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} 00 ${chrom_ginv_base} ${win_ginv_base} ${pheno_file} h0

        asreml -NS5 ${new_prefix}.as

        python ${baseDir}/bin/extract_loglik.py ${new_prefix}.asr none ${new_prefix}.llik ${chrom} ${window_size} asreml


        """
            
        }
}
