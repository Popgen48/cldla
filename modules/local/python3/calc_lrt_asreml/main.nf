process PYTHON3_CALC_LRT_ASREML {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir("${params.outdir}/python3/calc_lrt/${chrom}/", mode:'copy')
    maxForks 1

    input:
        tuple val(meta), path(chrom_giv), path(pheno_file), path(par_file), path(vcf_file)

    output:
        tuple val(meta), path("${outprefix}_${chrom}_results.txt"), emit: real_txt
        tuple val(meta), path("${outprefix}_${chrom}_perm_results.txt"), emit: perm_txt

    when:
        task.ext.when == null || task.ext.when

    script:
        chrom = meta.id
        outprefix = params.output_prefix
        window_size = params.window_size
        tool = params.tool
        n_perm = params.n_perm

        """
    python3 ${baseDir}/bin/vcf_to_local_lrt.py -v ${vcf_file} -r ${chrom} -w ${window_size} -c ${task.cpus} -g ${chrom_giv} -p ${pheno_file} -a ${par_file} -t ${tool} -o ${outprefix} -n ${n_perm}

        """
}
