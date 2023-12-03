process RUN_AIREMLF90_SERIAL{

    tag { "running_airemlf90_${chrom}" }
    label "process_medium"
    //maxForks 1
    publishDir("${params.outdir}/airemlf90_serial/${a_val}/", mode:"copy")

    input:
        tuple val(chrom), file(hap), file(map_f), file(win_ginv), file(chrom_ginv), file(pheno_file)
        val(a_val)

    output:
        path ("*.zip" ), emit:chrom_zip
        path ("*.sorted.lrt.out"), emit:lrt_out
        
    
    script:
        
        echidna_params = params.echidna_params
        template_h0 = a_val=="permute" ? params.pheno_file : pheno_file
        window_size = params.window_size
        

        
        """
        chrom_array=()
        
        for gin in \$(ls *.00.giv)
            do
                outgin=\$(basename \$gin .giv)

                python ${baseDir}/bin/prepare_airemlf90_params.py ${echidna_params} na \${gin} na ${template_h0} h0

                ${baseDir}/bin/blupf90+ \${outgin}.params

                tmp_arr=(\${gin//./ })

                chrom_array+=(\${tmp_arr[-3]})

                awk '\$0~/logL/{match(\$0,/(-2logL)([^0-9]+)([0-9\\.]+)(.*)/,a);print a[3]}' blupf90.log > \${outgin}.h0.llik

                cp blupf90.log \${outgin}.blupf90.log

                cp solutions \${outgin}.solutions

                #python ${baseDir}/bin/extract_loglik.py \${outgin}.asr none \${outgin}.h0.llik ${window_size}

            done                

        for h in \$(ls *.Hap)
            do

                diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' \${h})

                cnt=0

                outprefix=\$(basename \$h .Hap)

                for ph in \$(ls *.phe)
                    do
                        if [[ "\${cnt}" -eq 0 ]]
                            then
                                outprefix_pheno=\$outprefix
                        else
                            outprefix_pheno=\${outprefix}"_"\${cnt}
                        fi

                        cnt=\$((cnt+1))

                        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' \${h} \${ph} > \${outprefix_pheno}.dat

                        ar=(\${h//./ })

                        python ${baseDir}/bin/prepare_airemlf90_params.py ${echidna_params} \${diplo} \${ar[0]}.\${ar[1]}.00.giv \${outprefix}.giv \${outprefix_pheno}.dat h1

                        ${baseDir}/bin/blupf90+ \${outprefix_pheno}.params
                
                        awk '\$0~/logL/{match(\$0,/(-2logL)([^0-9]+)([0-9\\.]+)(.*)/,a);print a[3]}' blupf90.log > \${outprefix_pheno}.h1.llik
                        
                        cp blupf90.log \${outprefix_pheno}.blupf90.log

                        cp solutions \${outprefix_pheno}.solutions

                        #python ${baseDir}/bin/extract_loglik.py \${outprefix}.asr \${outprefix}.Map \${outprefix_pheno}.h1.llik ${window_size}

                    done

                done    
        
        for chrom in \${chrom_array}
            do
                mkdir \${chrom}

                mv *.\${chrom}.*.{log,solutions,llik} ./\${chrom}/

                cd ./\${chrom}

                cat *\${chrom}.*h1.llik > \${chrom}_window_h1.llik
            
                #python ${baseDir}/bin/calc_lrt.py \${chrom} *.h0.llik *h1.llik

                #cat *.lrt.out|awk '{printf "%d\\t%.0f\\t%.3f\\n",\$1,\$2,\$3}'|sort -n -k2,2 > \${chrom}.sorted.lrt.out

                #cp \${chrom}.sorted.lrt.out ../
    
                cd ..

                zip -r \${chrom}.airemlf90.out.zip ./\${chrom}/

            done


        """ 
}
