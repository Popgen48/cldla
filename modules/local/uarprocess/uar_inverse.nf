process UAR_INVERSE{

    tag { "uar_inverse_${chrom}" }
    label "oneCpu"
    //container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/uar_matrix/inverse/${chrom}", mode:"copy")


    input:
        tuple val(chrom), path(lgen), path(uar)

    output:
        tuple val(chrom), path ("*.giv"), emit: chrom_ginv
        
    
    script:

        dataset_id = params.output_prefix
        
        if( !params.run_wombat ){
        
        """
        ${baseDir}/bin/Bend5 ${uar} ${dataset_id}.${chrom}.B.grm

        python3 ${baseDir}/bin/ginverse.py ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.00.giv

        #python3 ${baseDir}/bin/bend_and_inverse.py ${uar} ${dataset_id}.${chrom}.00.giv
        #${baseDir}/bin/bend_and_inverse ${uar} ${dataset_id}.${chrom}.00.giv
        #Rscript ${baseDir}/bin/bend.R ${uar} ${dataset_id}.${chrom}.B.grm
        #python3 ${baseDir}/bin/ginverse.py ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.giv
        #diplo=\$(awk -v sc=0 '{if(!(\$1 in pop_id)){pop_id[\$1];sc++}}END{print sc}' ${lgen})
        #ginverse \${diplo} ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.giv 


        """ 
        }
        else{

        """
        ${baseDir}/bin/bend_and_inverse ${uar} ${dataset_id}.${chrom}.giv
        #Rscript ${baseDir}/bin/bend.R ${uar} ${dataset_id}.${chrom}.B.grm
        #python3 ${baseDir}/bin/ginverse.py ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.giv
        
        #Bend5 ${uar} ${dataset_id}.${chrom}.B.grm
        #diplo=\$(awk -v sc=0 '{if(!(\$1 in pop_id)){pop_id[\$1];sc++}}END{print sc}' ${lgen})
        #ginverse \${diplo} ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.ginv 
        cp .command.log ${dataset_id}.${chrom}.log
        python ${baseDir}/bin/prepare_ginv_wombat.py ${dataset_id}.${chrom}.log ${dataset_id}.${chrom}.ginv ${dataset_id}.${chrom}.giv


        """ 

        }
}
