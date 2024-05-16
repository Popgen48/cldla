process PYTHON3_CREATE_GCTA_INPUT{

    tag { "${tool}" }
    label "process_single"
    conda 'conda-forge::python=3.10'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://popgen48/python_bash:3.10-alpine' :
        'popgen48/python_bash:3.10-alpine' }"
    publishDir("${params.outdir}/python3/create_gcta_input/", mode:"copy")

    input:
        path(pheno)

    output:
        path ("*.phe"), emit: gcta_phe
        path ("*.qcovar"), emit: gcta_qcovar, optional: true
        path ("*.covar"), emit: gcta_covar, optional: true
        
    
    script:
        output_prefix = params.output_prefix

        """

        python3 ${baseDir}/bin/prepare_gcta_input.py ${pheno} ${output_prefix}

        """ 
}
