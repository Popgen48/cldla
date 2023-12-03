process CALC_LRT_RATIO{

    tag { "calculating_lrt_ratio_${chrom}" }
    label "process_single"
    container "popgen48/cldla_rpackages:1.0.0"
    publishDir("${params.outdir}/${tool}/LRT_ratio/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(h0_lrt), path(h1_lrt)

    output:
        path ( "*lrt.txt" ), emit: lrt
        
    
    script:
        //outprefix = if(params.echidna || params.asreml) ? h1_lrt.getName():
        tool= (params.echidna || params.asreml) ? (params.echidna ? "echidna":"asreml") : "airemlf90_serial"

        if(params.echidna || params.asreml ){
        
        """
        awk 'NR==FNR{nh[NR]=\$3;next}{if(\$3!="na"){print \$1,\$2,-2*(nh[FNR]-\$3)}}' ${h0_lrt} ${h1_lrt} > ${outprefix}.lrt

        """ 
        }
        else{
        
        """
        for z in \$(ls *map*h1.llik);do awk 'NR==FNR{h0=\$4;next}{print \$1,int(\$3),(h0-\$4)<0 ? 0:h0-\$4}' *00_map.h0.llik \$z;done|sort -n -k2,2 > chr${chrom}_cord_lrt.txt

        """

        }
}
