process RANDOMPHENO_CREATE{

    tag { "creating_random_phenodata_${chrom}" }
    label "oneCpu"
    publishDir("${params.outdir}/permutation_test/create_randomized_datasets/pheno/${chrom}", mode:"copy")

    input:
        tuple val(chrom), path(pheno)

    output:
        path ("*.phe"), emit: pheno_r
        
    
    script:
        
        def pheno_col = params.col_pheno
        def n_random_window = params.n_random_window


        """
        python3 ${baseDir}/bin/generate_pseudo_pheno_files.py ${pheno} ${chrom} ${pheno_col} ${n_random_window}
        
        """ 
}
