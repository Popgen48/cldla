process CREATE_HAPMAP {
    tag { "create hapmap_${chrom}_${outprefix}" }
    label 'process_single'
    container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/createhaplo/hapmap/${chrom}", mode:'copy')

    input:
        tuple val(chrom), val(window), val(outprefix), path(vcf)

    output:
        tuple val("${chrom}"), path('*.Hap'), emit: chrom_hap
        tuple val("${chrom}"), path('*.Map'), emit: chrom_map
        tuple val("${chrom}"), path('*.par'), emit: chrom_par

    script:

        window_size = params.window_size
        dataset_id = params.output_prefix

        """

        python3 ${baseDir}/bin/vcf_to_hapmap.py ${vcf} ${chrom} ${window_size} ${window} ${dataset_id} ${outprefix}

        """
}
