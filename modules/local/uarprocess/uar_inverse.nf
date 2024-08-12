process UAR_INVERSE {
    tag { "uar_inverse_${chrom}" }
    label 'process_single'
    //container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/uar_matrix/inverse/${chrom}", mode:'copy')

    input:
        tuple val(chrom), path(lgen), path(uar)

    output:
        tuple val(meta), path('*.giv'), emit: giv

    script:

        dataset_id = params.output_prefix
        meta = [id:chrom]

        """
        ${baseDir}/bin/bend ${uar} ${dataset_id}.${chrom}.B.grm

        diplo=\$(awk -v sc=0 '{if(!(\$1 in pop_id)){pop_id[\$1];sc++}}END{print sc}' ${lgen})

        ${baseDir}/bin/ginverse \${diplo} ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.giv

        """
}
