process PYTHON3_FILTER_VCF {
    tag { "${meta.id}" }
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py310h41dec4a_0' :
        'quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0' }"
    publishDir("${params.outdir}/python3/filter_vcf/", mode:'copy')

    input:
        tuple val(meta), path(vcf), path(idx), path(pheno_file)

    output:
        tuple val(meta), path("${outprefix}.vcf.gz"), emit: gzvcf

    when:
        task.ext.when == null || task.ext.when

    script:
        chrom = meta.id
        outprefix = params.output_prefix + '.' + chrom
        maf = params.maf

        """

        python3 ${baseDir}/bin/filter_vcf.py -v ${vcf} -s ${pheno_file} -m ${maf} -o ${outprefix}.vcf.gz

        """
}
