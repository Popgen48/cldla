process PYTHON3_PREPARE_PARAMETER_TEMPLATE{

    tag { "${tool}" }
    label "process_single"
    conda 'conda-forge::python=3.10'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://popgen48/python_bash:3.10-alpine' :
        'popgen48/python_bash:3.10-alpine' }"
    publishDir("${params.outdir}/python3/prepare_parameter_template/", mode:"copy")

    input:
        path (pheno)
        path (add_params)

    output:
        path ("*.params"), emit: blp_template, optional: true
        path ("*.as"), emit: asr_template, optional: true
        
    
    script:
        output_prefix = params.output_prefix
        tool = params.tool
        
        

        """

        python3 ${baseDir}/bin/prepare_parameter_templates.py ${pheno} ${tool} ${add_params} ${output_prefix}

        """ 
}
