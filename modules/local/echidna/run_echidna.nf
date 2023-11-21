process RUN_ECHIDNA{

    tag { "running_echidna_${chrom}" }
    label "process_medium"
    maxForks 1
    publishDir("${params.outdir}/echidna/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(hap), path(map_f), path(win_ginv), path(chrom_ginv), path(pheno_file)

    output:
        tuple val("${chrom}"), path ("*.{as,esr}"), emit: chrom_as
        tuple val("${chrom}"), path ( "*.llik" ), emit: chrom_llik
        
    
    script:
        
        echidna_params = params.echidna_params
        win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        chrom_ginv_base = chrom_ginv.getName()
        new_prefix = hap==[] ? chrom_ginv.getName().minus('.giv'):hap.getName().minus('.Hap')

        if( hap != [] ){
        
        """

        diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' ${hap})

        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' ${hap} ${pheno_file} > ${new_prefix}.phe
        
        python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} \${diplo} ${chrom_ginv_base} ${win_ginv_base} ${new_prefix}.phe h1

        ${baseDir}/bin/Echidna -w2ro ${new_prefix}.as

        cp ./${new_prefix}/*.esr .

        python ${baseDir}/bin/extract_loglik.py ${new_prefix}.esr ${map_f} ${new_prefix}.llik


        """ 
        }
        else{

        """

        python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} 00 ${chrom_ginv_base} ${win_ginv_base} ${pheno_file} h0

        ${baseDir}/bin/Echidna -w2ro ${new_prefix}.as
        
        cp ./${new_prefix}/*.esr .

        python ${baseDir}/bin/extract_loglik.py ${new_prefix}.esr none ${new_prefix}.llik

        """
            
        }
}
