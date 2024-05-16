process H2_RANDOMPHENO_CREATE{

    tag { "creating_random_phenodata" }
    label "process_single"
    container "popgen48/cldla_python_r_packages:1.0.0"
    publishDir("${params.outdir}/gcta/ranomized_pheno_datasets/", mode:"copy")

    input:
        path(pheno)

    output:
        path ("*.phe"), emit: pheno_he
        
    
    script:
        
        def n_random_window = params.n_perm


        """
        python3 ${baseDir}/bin/generate_pseudo_pheno_h2.py ${pheno} ${n_random_window}
        
        """ 
}
