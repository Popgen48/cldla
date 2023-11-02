process RUN_ECHIDNA{

    tag { "running_echidna_${chrom}" }
    label "oneCpu"
    //container "popgen48/cldla_rpackages:1.0.0"
    //conda "${baseDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
    //    'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/echidna/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(hap), path(win_ginv), path(chrom_ginv), path(pheno_file)

    output:
        tuple val("${chrom}"), path ("*.{as,esr}"), emit: chrom_as
        
    
    script:
        
        echidna_params = params.echidna_params
        win_ginv_base = win_ginv.getName()
        chrom_ginv_base = chrom_ginv.getName()
        new_prefix = hap.getName().minus('.Hap')
        
        """
        diplo=\$(awk 'END{print \$4}' ${hap})

        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' ${hap} ${pheno_file} > ${new_prefix}.phe
        
        python ${baseDir}/bin/prepare_echidna_params.py ${echidna_params} \${diplo} ${chrom_ginv_base} ${win_ginv_base} ${new_prefix}.phe

        ${baseDir}/bin/Echidna -w2ro ${new_prefix}.as

        """ 
}
