process UAR_INVERSE{

    tag { "uar_inverse_${chrom}" }
    label "process_single"
    container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/uar_matrix/inverse/${chrom}", mode:"copy")


    input:
        tuple val(chrom), path(lgen), path(uar)

    output:
        tuple val(chrom), path ("*.giv"), emit: chrom_ginv
        
    
    script:

        dataset_id = params.output_prefix
        
        
        """
        ${baseDir}/bin/Bend5 ${uar} ${dataset_id}.${chrom}.B.grm

        python3 ${baseDir}/bin/ginverse.py ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.00.giv



        """ 
}
