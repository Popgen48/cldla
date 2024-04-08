process PLOT_LRT_RATIO{

    tag { "plot_lrt_ratio" }
    label "process_single"
    container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/${tool}/", mode:"copy")

    input:
        path(lrt_files)

    output:
        path ( "*.html" )
        path ("genomewide_lrt.out")
        
    
    script:
        outprefix = params.output_prefix
        lrt_threshold = params.lrt_threshold
        tool= (params.echidna || params.asreml) ? (params.echidna ? "echidna":"asreml") : "airemlf90_serial"
        
        """
        cat ${lrt_files} | sort -n -k1,1 -k2,2 > genomewide_lrt.out

        python3 ${baseDir}/bin/make_manhattan_plot.py genomewide_lrt.out ${baseDir}/manhattanplot.yml ${lrt_threshold} 100 LRT_values ${outprefix}

        """ 
}
