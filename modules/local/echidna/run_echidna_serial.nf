process RUN_ECHIDNA_SERIAL{

    tag { "running_echidna_serial" }
    label "process_medium"
    maxForks 1
    publishDir("${params.outdir}/echidna/", mode:"copy")

    input:
        tuple val(chrom), file(hap), file(map_f), file(win_ginv), file(chrom_ginv), file(pheno_file)

    output:
        //path ("*.{as,esr}"), emit: chrom_as
        path ("*.zip" ), emit:chrom_zip
        path ("*.lrt.out"), emit:lrt_out
        //path ( "*.llik" ), emit: chrom_llik
        
    
    script:
        
        echidna_params = params.echidna_params
        //pheno_file = params.pheno_file
        //win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        //new_prefix = hap==[] ? chrom+"_00":hap.getName().minus('.Hap')

        
        """
        chrom_array=()
        
        for gin in \$(ls *.00.giv)
            do
                outgin=\$(basename \$gin .giv)

                python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} na \${gin} na ${pheno_file} h0

                ${baseDir}/bin/Echidna -w2ro \${outgin}.as

                tmp_arr=(\${gin//./ })

                chrom_array+=(\${tmp_arr[-3]})

                cat \${gin} > ./\${outgin}/\${outgin}.giv

                python ${baseDir}/bin/extract_loglik.py ./\${outgin}/\${outgin}.esr none \${outgin}.h0.llik

                zip -r \${outgin}.zip ./\${outgin}/

                rm *.{esk,eqo}

                rm -r ./\${outgin}

            done                

        for h in \$(ls *.Hap)
            do
                outprefix=\$(basename \$h .Hap)

                diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' \${h})

                awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' \${h} ${pheno_file} > \${outprefix}.phe

                ar=(\${h//./ })

                python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} \${diplo} \${ar[0]}.\${ar[1]}.00.giv \${outprefix}.giv \${outprefix}.phe h1

                ${baseDir}/bin/Echidna -w2ro \${outprefix}.as
                
                rm *.{esk,eqo}

                cat \${outprefix}.giv > ./\${outprefix}/\${outprefix}.giv

                cat \${outprefix}.phe > ./\${outprefix}/\${outprefix}.phe

                python ${baseDir}/bin/extract_loglik.py ./\${outprefix}/\${outprefix}.esr \${outprefix}.Map \${outprefix}.h1.llik

                zip -r \${outprefix}.zip ./\${outprefix}/

                rm -r ./\${outprefix}/

            done    
        
        for chrom in \${chrom_array}
            do

                mkdir \${chrom}

                mv *.\${chrom}.*.zip ./\${chrom}/

                mv *\${chrom}.*h1.llik ./\${chrom}/

                mv *\${chrom}.*h0.llik ./\${chrom}/

                #mv *.\${chrom}.*.{llik,as,asr,res,sln} ./\${chrom}/

                cd ./\${chrom}

                cat *\${chrom}.*h1.llik > \${chrom}_window_h1.llik
            
                python ${baseDir}/bin/calc_lrt.py \${chrom} *.h0.llik *h1.llik

                cp *.lrt.out ../
    
                cd ..

                zip -r \${chrom}.echidna.out.zip ./\${chrom}/

                rm -r ./\${chrom}/

            done

        #unset chrom_array

        #unset tmp_arr

        """ 
}
