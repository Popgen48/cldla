process PYTHON3_PREPARE_MANHATTAN_PLOT_INPUT {
    tag { "${outprefix}" }
    label 'process_medium'
    publishDir("${params.outdir}/python3/prepare_manhattan_plot_input/", mode:'copy')
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://popgen48/prepare_manhattan_input:1.0.0' :
        'popgen48/prepare_manhattan_input:1.0.0' }"

    input:
        path(chrom_vcf_idx)
        path(lrt_txt)

    output:
        path("*_maninp.txt"), emit: maninp_txt optional true
        path("*_threshold.txt"), emit: threshold_txt optional true

    when:
        task.ext.when == null || task.ext.when

    script:
        outprefix = params.output_prefix
        def args = task.ext.args ?: ''

        """

        python3 ${baseDir}/bin/prepare_manhattan_plot_inputs.py ${chrom_vcf_idx} ${outprefix} ${args} ${lrt_txt}

        """
}
