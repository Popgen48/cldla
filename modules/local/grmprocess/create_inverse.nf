process CREATE_INVERSE{

    tag { "creating grm_and_inverse" }
    label "oneCpu"
    //conda "${baseDir}/environment.yml"
    //container 'popgen48/cldla_python_r_packages:1.0.0'
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
        if ( !params.run_wombat ){

        """
        
        ${baseDir}/bin/cLDLA_snp ${o_prefix}
        
        ${baseDir}/bin/Bend5 ${o_prefix}.grm ${o_prefix}.B.grm

        python3 ${baseDir}/bin/ginverse.py ${o_prefix}.B.grm ${o_prefix}.giv

        #python3 ${baseDir}/bin/bend_and_inverse.py ${o_prefix}.grm ${o_prefix}.giv

        #${baseDir}/bin/bend_and_inverse ${o_prefix}.grm ${o_prefix}.giv

        #Rscript ${baseDir}/bin/bend.R ${o_prefix}.grm ${o_prefix}.B.grm


        #Bend5 ${o_prefix}.grm ${o_prefix}B.grm
        #diplo=\$(awk -v sc=0 '{if(!(\$4 in pop_id)){pop_id[\$4];sc++}}END{print sc}' ${hap})
        #ginverse \${diplo} ${o_prefix}B.grm ${o_prefix}.giv

        """ 
        }
    else{
        
        """
        
        cLDLA_snp ${o_prefix}

        #Rscript ${baseDir}/bin/bend.R ${o_prefix}.grm ${o_prefix}.B.grm

        #python3 ${baseDir}/bin/ginverse.py ${o_prefix}.B.grm ${o_prefix}.giv

        #Bend5 ${o_prefix}.grm ${o_prefix}B.grm
        #diplo=\$(awk -v sc=0 '{if(!(\$4 in pop_id)){pop_id[\$4];sc++}}END{print sc}' ${hap})
        #ginverse \${diplo} ${o_prefix}B.grm ${o_prefix}.ginv
        #cp .command.log ${o_prefix}.log
        python ${baseDir}/bin/prepare_ginv_wombat.py ${o_prefix}.log ${o_prefix}.ginv ${o_prefix}.giv


        """ 
    }
}
