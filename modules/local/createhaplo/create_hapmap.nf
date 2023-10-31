process CREATE_HAPMAP{

    tag { "create hapmap_${chrom}_${window}" }
    label "oneCpu"
    //conda "${baseDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
        'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/createhaplo/hapmap/${chrom}", mode:"copy")

    input:
        tuple val(chrom), val(window), path(vcf)

    output:
        tuple val("${chrom}"), path ("*.Hap"), emit: chrom_hap
        tuple val("${chrom}"), path ("*.Map"), emit: chrom_map
        tuple val("${chrom}"), path ("*.par"), emit: chrom_par
        
        
    
    script:
        
        window_size = params.window_size
        dataset_id = params.output_prefix

        """
        
        python3 ${baseDir}/bin/vcf_to_hapmap.py ${vcf} ${chrom} ${window_size} ${window} ${dataset_id}

        """ 
}
