process MAKE_GRM{

    tag { "make grm" }
    label "process_medium"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/gcta/make_grm/", mode:"copy")

    input:
        tuple val(chrom), path(bim)
        tuple val(chrom), path(bed)
        tuple val(chrom), path(fam)

    output:
        path ("*.grm.id"), emit: gcta_grm_id
        path ("*.grm.N.bin"), emit: gcta_grm_n_bin
        path ("*.grm.bin"), emit: gcta_grm_bin
        
    
    script:
    
        i_prefix = bed.getBaseName().minus('.bed')
        num_autosomes = params.num_autosomes

        """
        gcta --bfile ${i_prefix} --autosome --make-grm --autosome-num ${num_autosomes} --out ${i_prefix} --thread-num ${task.cpus}

        """ 
}
