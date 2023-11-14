process CALC_LRT_RATIO{

    tag { "calculating_lrt_ratio_${chrom}" }
    label "oneCpu"
    //container "popgen48/cldla_rpackages:1.0.0"
    //conda "${baseDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
    //    'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/echidna/LRT_ratio/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(h0_lrt), path(h1_lrt)

    output:
        tuple val("${chrom}"), path ( "*.lrt" ), emit: chrom_lrt
        
    
    script:
        outprefix = h1_lrt.getName()
        
        """
        awk 'NR==FNR{nh[NR]=\$3;next}{if(\$3!="na"){print \$1,\$2,-2*(nh[FNR]-\$3)}}' ${h0_lrt} ${h1_lrt} > ${outprefix}.lrt

        """ 
}
