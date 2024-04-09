process PYTHON3_VCFTOLGEN{

    tag { "creating_lgen_${chrom}" }
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py310h41dec4a_0':
        'quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0' }"
    publishDir("${params.outdir}/python3/filter_vcf/", mode:"copy")

    input:
        tuple val(chrom), path(vcf), path(pheno_file)

    output:
        tuple val("${chrom}"), path ("*.lgen"), emit: chrom_lgen
        
    
    script:
        
        dataset_id = params.output_prefix

        """
        
        python3 ${baseDir}/bin/vcf_to_lgen.py ${vcf} ${chrom} ${pheno_file} ${dataset_id}_no_${chrom}

        """ 
}
