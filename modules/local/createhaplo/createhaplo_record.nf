process CREATEHAPLO_RECORD{

    tag { "create haplo_record_${chrom}" }
    label "process_single"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/createhaplo/record/${chrom}/", mode:"copy")

    input:
        tuple val(chrom), path(vcf)

    output:
        path ("*window_counts*"), emit: record
        path ("*random_windows*"), emit: random_record
        
    
    script:
        def window_size = params.window_size
        def n_random_window = params.n_random_window

        """
        
        python3 ${baseDir}/bin/get_vcf_record_count.py ${vcf} ${chrom} ${window_size} ${n_random_window}

        """ 
}
