process CALC_LRT_RATIO{

    tag { "calculating_lrt_ratio_${chrom}" }
    label "process_single"
    container "popgen48/cldla_rpackages:1.0.0"
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
