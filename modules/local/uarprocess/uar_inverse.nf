process UAR_INVERSE{

    tag { "uar_inverse_${chrom}" }
    label "oneCpu"
    container 'popgen48/cldla_binaries:1.0.0'
    publishDir("${params.outdir}/uar_matrix/inverse/${chrom}", mode:"copy")


    input:
        tuple val(chrom), path(lgen), path(uar)

    output:
        tuple val(chrom), path ("*.giv"), emit: chrom_ginv
        
    
    script:

        dataset_id = params.output_prefix


        """
        
        Bend5 ${uar} ${dataset_id}.${chrom}.B.grm
        diplo=\$(awk -v sc=0 '{if(!(\$1 in pop_id)){pop_id[\$1];sc++}}END{print sc}' ${lgen})
        ginverse \${diplo} ${dataset_id}.${chrom}.B.grm ${dataset_id}.${chrom}.giv 


        """ 
}
