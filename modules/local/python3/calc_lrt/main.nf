process PYTHON3_CALC_LRT{

    tag { "${meta.id}" }
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py310h41dec4a_0':
    //    'quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0' }"
    publishDir("${params.outdir}/python3/calc_lrt/${chrom}/", mode:"copy")
    maxForks 1

    input:
        tuple val(meta), path(chrom_giv), path(pheno_file), path(par_file), path(vcf_file)

    output:
        tuple val(meta), path("*_results.txt"), emit: txt

    when:
        task.ext.when == null || task.ext.when

    script:
        chrom = meta.id
        outprefix = params.output_prefix
        window_size = params.window_size
        num_process = params.n_var_est_process
        tool = params.tool
        n_perm = params.n_perm

        """
    python3 ${baseDir}/bin/vcf_to_local_lrt.py -v ${vcf_file} -r ${chrom} -w ${window_size} -c ${task.cpus} -g ${chrom_giv} -p ${pheno_file} -a ${par_file} -t ${tool} -o ${outprefix} -n ${n_perm}


        """ 

}
