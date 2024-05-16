process REML_GRM{

    tag { "reml grm" }
    label "process_medium"
    conda "bioconda::gcta==1.94.1"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/gcta/reml_grm/", mode:"copy")

    input:
        tuple path(grm_i), path(grm_nb), path(grm_b), path(pheno_f), path(qcovar_f), path(covar_f)

    output:
        path ("*.hsq"), emit: gcta_hsq
        
    
    script:
        def args1 = qcovar_f!=[] ? " --qcovar "+qcovar_f : ""
        def args2 = covar_f!=[] ? " --covar "+covar_f : ""
        output_prefix = pheno_f.getBaseName().minus(".phe")
        grm_f = grm_i.getName().minus(".grm.id")


        """
        gcta --reml --grm ${grm_f} --pheno ${pheno_f} --out ${output_prefix} ${args1} ${args2} --thread-num ${task.cpus}

        """ 
}
