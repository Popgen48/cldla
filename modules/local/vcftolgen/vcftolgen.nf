process VCFTOLGEN{

    tag { "creating_lgen_${chrom}" }
    label "oneCpu"
    //conda "${baseDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
        'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/vcftolgen/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(vcf)

    output:
        tuple val("${chrom}"), path ("*.lgen"), emit: lgen
        
    
    script:
        
        dataset_id = params.output_prefix

        """
        
        python3 ${baseDir}/bin/vcf_to_lgen.py ${vcf} ${chrom} ${dataset_id}_no_${chrom}

        """ 
}
