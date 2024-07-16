process PYTHON3_MANHATTAN_PLOT{

    tag { "${outprefix}" }
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://popgen48/plot_selection_results:1.0.0' :
        'popgen48/plot_selection_results:1.0.0' }"
    publishDir("${params.outdir}/python3/manhattan_plot/", mode:"copy")

    input:
        path(merged_result)
        path(yml)
        path(chrom_threshold)

    output:
        path("*.html"), emit: html

    script:
        outprefix = params.output_prefix
        def args = task.ext.args ?: ''

        """
        
        python ${baseDir}/bin/make_manhattan_plot.py ${merged_result} ${yml} ${chrom_threshold} ${args} ${outprefix}
        

        """ 
}
