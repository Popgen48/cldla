process RUN_WOMBAT{

    tag { "running_wombat_${chrom}" }
    label "process_single"
    //container "popgen48/cldla_rpackages:1.0.0"
    //conda "${baseDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
    //    'biocontainers/pysam:v0.15.2ds-2-deb-py3_cv1' }"
    publishDir("${params.outdir}/wombat/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(hap), path(map_f), path(win_ginv), path(chrom_ginv), path(pheno_file)

    output:
        tuple val("${chrom}"), path ("Iterates"), emit: chrom_as
        tuple val("${chrom}"), path ( "WOMBAT.log" ), emit: chrom_llik
        
    
    script:
        
        wombat_params = params.wombat_params
        win_ginv_base = hap==[] ? chrom+"_00.00": win_ginv.getName()
        chrom_ginv_base = chrom_ginv.getName()
        new_prefix = hap==[] ? chrom+"_00":hap.getName().minus('.Hap')

        if( hap != [] ){
        
        """

        diplo=\$(awk -v max=0 '{if(\$4>max){max=\$4;next}else;next}END{print max}' ${hap})

        awk 'NR==FNR{sample[\$1]=\$4;next}{\$(NF+1)=sample[\$1];print}' ${hap} ${pheno_file} > ${new_prefix}.phe

        cat ${win_ginv} > Ibd.gin

        cat ${chrom_ginv} > Indi.gin
        
        python ${baseDir}/bin/prepare_wombat_params.py ${wombat_params} \${diplo} ${new_prefix}.phe h1

        ${baseDir}/bin/wombat wombat.par 

        #cp ./${new_prefix}/*.esr .

        #python ${baseDir}/bin/extract_loglik.py ${new_prefix}.esr ${map_f} ${new_prefix}.llik


        """ 
        }
        else{

        """
        cat ${chrom_ginv} > Indi.gin

        python ${baseDir}/bin/prepare_wombat_params.py ${wombat_params} \${diplo} ${new_prefix}.phe h0

        ${baseDir}/bin/wombat wombat.par
        
        #cp ./${new_prefix}/*.esr .

        #python ${baseDir}/bin/extract_loglik.py ${new_prefix}.esr none ${new_prefix}.llik

        """
            
        }
}
