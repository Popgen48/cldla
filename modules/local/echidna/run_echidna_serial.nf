process RUN_ECHIDNA_SERIAL{

    tag { "running_echidna_serial" }
    label "process_medium"
    //container "popgen48/cldla_rpackages:1.0.0"
    //conda "${baseDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
    //    'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/asreml_serial/", mode:"copy")

    input:
        path(hap)
        path(map_f)
        path(win_ginv)
        path(chrom_ginv)
        path(pheno_file)

    output:
        //path ("*.{as,esr}"), emit: chrom_as
        path ("*.zip" ), emit:chrom_zip
        path ("*.lrt.out"), emit:lrt_out
        //path ( "*.llik" ), emit: chrom_llik
        
    
    script:
        
        echidna_params = params.echidna_params
        //win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        //new_prefix = hap==[] ? chrom+"_00":hap.getName().minus('.Hap')

        if( hap != [] ){
        
        """
        chrom_array=()
        
        for gin in \$(ls *.00.giv)
            do
                outgin=\$(basename \$gin .giv)

                python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} na \${gin} na ${pheno_file} h0

                asreml -NS5 \${outgin}.as

                tmp_arr=(\${gin//./ })

                chrom_array+=(\${tmp_arr[-3]})

                python ${baseDir}/bin/extract_loglik.py \${outgin}.asr none \${outgin}.h0.llik

            done                

        for h in \$(ls *.Hap)
            do
                outprefix=\$(basename \$h .Hap)

                diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' \${h})

                awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' \${h} ${pheno_file} > \${outprefix}.phe

                ar=(\${h//./ })

                python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} \${diplo} \${ar[0]}.\${ar[1]}.00.giv \${outprefix}.giv \${outprefix}.phe h1

                asreml -NS5 \${outprefix}.as
                
                rm *.{tsv,tmp,ask,veo,rsv,msv,vvp,yht}

                python ${baseDir}/bin/extract_loglik.py \${outprefix}.asr \${outprefix}.Map \${outprefix}.h1.llik

            done    
        
        for chrom in \${chrom_array}
            do
                mkdir \${chrom}

                mv *.\${chrom}.*.{llik,as,asr,res,sln} ./\${chrom}/

                cd ./\${chrom}

                cat *\${chrom}.*h1.llik > \${chrom}_window_h1.llik
            
                python ${baseDir}/bin/calc_lrt.py \${chrom} *.h0.llik *h1.llik

                cp *.lrt.out ../
    
                cd ..

                zip -r \${chrom}.asreml.out.zip ./\${chrom}/

            done

        #unset chrom_array

        #unset tmp_arr

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
