process PYTHON3_CALC_LRT_BLUPF90 {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir("${params.outdir}/python3/calc_lrt/${chrom}/", mode:'copy')
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://popgen48/cldla:1.0.0' :
        'popgen48/cldla:1.0.0'}"
    maxForks 1

    input:
        tuple val(meta), path(chrom_giv), path(pheno_file), path(par_file), path(vcf_file)

    output:
        tuple val(meta), path("*filtered_window_results.txt"), emit: real_txt
        tuple val(meta), path("*all_window_results.txt"), emit: all_window
        tuple val(meta), path("*perm_results.txt"), emit: perm_txt
        path("*.csv") optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        chrom = meta.id
        outprefix = params.output_prefix
        window_size = params.window_size
        tool = params.tool
        n_perm = params.n_perm

        if(params.store){
                args = args+ " --s "+" -O "+ params.store_outdir
            }
        """
    python3 ${baseDir}/bin/vcf_to_local_lrt.py -v ${vcf_file} -r ${chrom} -w ${window_size} -c ${task.cpus} -g ${chrom_giv} -p ${pheno_file} -a ${par_file} -t ${tool} -o ${outprefix} -n ${n_perm} ${args}

        """
}
