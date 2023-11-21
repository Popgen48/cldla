process RUN_ASREML_SERIAL{

    tag { "running_asreml_serial" }
    label "process_medium"
    maxForks 1
    publishDir("${params.outdir}/asreml_serial/", mode:"copy")

    input:
        tuple val(chrom), file(hap), file(map_f), file(win_ginv), file(chrom_ginv)

    output:
        path ("*.zip" ), emit:chrom_zip
        path ("*.lrt.out"), emit:lrt_out
        
    
    script:
        
        echidna_params = params.echidna_params
        pheno_file = params.pheno_file

        
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


        """ 
}
