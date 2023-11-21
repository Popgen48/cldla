process CALC_UAR{

    tag { "calculating_uar_${chrom}" }
    label "process_single"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/uar_matrix/create/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(lgen)

    output:
        tuple val("${chrom}"), path ("*.grm"), emit: chrom_uar
        
    
    script:
        
        dataset_id = params.output_prefix

        """
        
        Rscript ${baseDir}/bin/calc_uar_matrix.r ${lgen} ${dataset_id}_${chrom}.uar.txt

        awk -v sum=0 'NR>1{sum++;for(i=2;i<=NF;i++){if(i>=NR){print sum,i-1,\$i}}}' ${dataset_id}_${chrom}.uar.txt > ${dataset_id}_${chrom}.pregrm.txt

        diplo=\$(awk -v sc=0 '{if(!(\$1 in pop_id)){pop_id[\$1];sc++}}END{print sc}' ${lgen})

        awk -v a=\$diplo 'BEGIN{print a}{print}' ${dataset_id}_${chrom}.pregrm.txt > ${dataset_id}_${chrom}.grm

        """ 
}
