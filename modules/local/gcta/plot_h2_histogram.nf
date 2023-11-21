process PLOT_H2_HISTOGRAM{

    tag { "reml grm" }
    label "process_medium"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/gcta/heritability_output/", mode:"copy")

    input:
        path(h2_ori)
        path(h2_sim)

    output:
        path ("*.{html,txt}")
        
    
    script:

        output_prefix = params.output_prefix

        """
        for z in \$(ls dataset*.hsq);do awk -v cnt=0 '{cnt++;if(cnt==5){print \$2}}' \${z} >> simulated_h2.txt;done

        awk -v cnt=0 '{cnt++;if(cnt==5){print \$2}}' ${h2_ori} > observed_h2.txt

        python3 ${baseDir}/bin/plot_histogram.py simulated_h2.txt observed_h2.txt ${output_prefix}

        """ 
}
