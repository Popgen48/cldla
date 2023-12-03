process RUN_AIREMLF90_WINDOW_WISE{

    tag { "running_airemlf90_${chrom}" }
    label "process_medium"
    //maxForks 1
    publishDir("${params.outdir}/airemlf90_serial/${chrom}/${a_val}/", mode:"copy")

    input:
        tuple val(chrom), file(hap), file(map_f), file(win_ginv), file(chrom_ginv), file(pheno_file)
        val(a_val)

    output:
        tuple val("${chrom}"), path ("*map*.llik"), emit: chrom_llik
        path ("*.solutions"), emit: solution
        
    
    script:
        
        window_size = params.window_size
        echidna_params = params.echidna_params
        win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        chrom_ginv_base = chrom_ginv.getName()
        new_prefix = hap==[] ? chrom_ginv.getName().minus('.giv'):hap.getName().minus('.Hap')
        template_h0 = a_val=="permute" ? params.pheno_file : pheno_file

        if( hap == [] ){
        
        """
        
        outgin=\$(basename ${chrom_ginv} .giv)

        python3 ${baseDir}/bin/prepare_airemlf90_params.py ${echidna_params} na ${chrom_ginv} na ${template_h0} h0

        ${baseDir}/bin/blupf90+ \${outgin}.params

        outprefix_pheno=\$(basename ${pheno_file} .phe)

        awk '\$0~/logL/{match(\$0,/(-2logL)([^0-9]+)([0-9\\.]+)(.*)/,a);print a[3]}' blupf90.log > \${outgin}.h0.llik

        python3 ${baseDir}/bin/extract_loglik.py \${outgin}.h0.llik none \${outgin}_map.h0.llik ${chrom} ${window_size} blupf90+

        cp blupf90.log \${outgin}.blupf90.log

        cp solutions \${outgin}.solutions

        """
        }

        else{
        
        """

        diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' ${hap})

        outprefix=\$(basename ${hap} .Hap)

        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' ${hap} ${pheno_file} > \${outprefix}.dat

        python ${baseDir}/bin/prepare_airemlf90_params.py ${echidna_params} \${diplo} ${chrom_ginv} ${win_ginv} \${outprefix}.dat h1

        ${baseDir}/bin/blupf90+ \${outprefix}.params
                
        awk '\$0~/logL/{match(\$0,/(-2logL)([^0-9]+)([0-9\\.]+)(.*)/,a);print a[3]}' blupf90.log > \${outprefix}.h1.llik

        python3 ${baseDir}/bin/extract_loglik.py \${outprefix}.h1.llik ${map_f} \${outprefix}_map.h1.llik ${chrom} ${window_size} blupf90+
                        
        cp blupf90.log \${outprefix}.blupf90.log

        cp solutions \${outprefix}.solutions


        """ 
    }
}
