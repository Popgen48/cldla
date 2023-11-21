process VCFTOLGEN{

    tag { "creating_lgen_${chrom}" }
    label "process_single"
    container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/vcftolgen/", mode:"copy")

    input:
        tuple val(chrom), path(vcf)

    output:
        tuple val("${chrom}"), path ("*.lgen"), emit: chrom_lgen
        
    
    script:
        
        dataset_id = params.output_prefix

        """
        
        python3 ${baseDir}/bin/vcf_to_lgen.py ${vcf} ${chrom} ${dataset_id}_no_${chrom}

        """ 
}
