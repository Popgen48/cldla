process CREATE_INVERSE{

    tag { "creating grm_and_inverse" }
    label "process_single"
    container 'popgen48/cldla_python_r_packages:1.0.0'
    publishDir("${params.outdir}/grm/inverse/", mode:"copy")

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
        
        ${baseDir}/bin/cLDLA_snp ${o_prefix}
        
        ${baseDir}/bin/Bend5 ${o_prefix}.grm ${o_prefix}.B.grm

        python3 ${baseDir}/bin/ginverse.py ${o_prefix}.B.grm ${o_prefix}.giv


        """ 
}
