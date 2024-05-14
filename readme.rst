Introduction
-------------

popgen-cldla is a fully automated workflow to identify QTLs using combined linkage disequilibrium and linkage analysis (CLDLA). 

| The workflow consists of the following steps:

| 1. Phase the genotypes using beagle or shapeit (optional)
| 2. Filtering the SNPs (optional)
| 3. Estimation of genetic relationship across genome (GRM)
| 4. Estimation of genetic relationship at each putative QTL region (DRM)
| 5. Bending and inverting GRM and DRM
| 6. Variance component estimation 
| 7. Plotting the results in the form of manhattan plot

Input files
-------------
| 1. List of phased or unphased vcf files in a csv format

..  code-block:: Bash
    
    10,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr10.TxT.vcf.gz,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr10.TxT.vcf.gz.csi
    11,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr11.TxT.vcf.gz,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr11.TxT.vcf.gz.csi
| 2. Phenotype file 
| 3. A template file of parameters with the definition of mixed and fixed effects.
