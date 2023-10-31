process CREATE_INVERSE{

    tag { "creating grm_and_inverse" }
    label "oneCpu"
    //conda "${baseDir}/environment.yml"
    container 'popgen48/cldla_binaries:1.0.0'
    publishDir("${params.outdir}/grm/inverse/${chrom}", mode:"copy")
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
    //    'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"

    input:
        tuple val(chrom), path(hap)
        tuple val(chrom), path(map)
        tuple val(chrom), path(par)

    output:
        tuple val(o_prefix), path ("*.giv"), emit: chrwin_ginv
        
    
    script:

        o_prefix = hap.getBaseName().minus('.Hap')
        //chrom = prefix.split(".")[0]


        """
        
        cLDLA_snp ${o_prefix}
        Bend5 ${o_prefix}.grm ${o_prefix}B.grm
        diplo=\$(awk 'END{print \$1}' ${hap})
        ginverse \${diplo} ${o_prefix}B.grm ${o_prefix}.giv
 


        """ 
}
