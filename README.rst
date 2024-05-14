Introduction
-------------

popgen-cldla is a fully automated workflow to identify QTLs using combined linkage disequilibrium and linkage analysis (CLDLA). Unlike genome-wide association study (GWAS), cLDLA is a mixed-linear-model (MLM) based apporach to estimate variance component of each putative QTL using the specialized tools such as ASREML and blupf90+. 

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
   .
   .
   .

| 2. Phenotype file 

..  code-block:: Bash
    
    1	MLS169	1	11	10.6	47.0	32.7
    2	MLS170	1	9	12.0	44.2	33.6
    .
    .
    .

| *Note that first, second and the last column in the phenotype file should be the index of individuals, individual id and the trait (Y), respectively.*

| 3. A template file of parameters with the definition of linear mixed models.

| Parameter template file of *asreml*
..  code-block:: Bash
	This line is a headline: TailMLS04 design include all fix effects and phenotype
	indi 362 !I
	TierLID 362 !A
	Sex 1
	Age 1
	BodyW 1
	WitherH 1
	TailL 1
	!I
	!LDET
	!LDET
	!AISING !NODISPLAY !MAXIT 99
	TailL ~ mu Sex Age BodyW WitherH !r giv(iDip,2) !r giv(indi,1)

| Parameter template file of *blupf90+*
..  code-block:: Bash
	NUMBER_OF_TRAITS
	1
	NUMBER_OF_EFFECTS
	5
	OBSERVATION(S)
	7
	WEIGHT(S)

	EFFECTS:
	1 362 cross
	3  2 cov
	4  41 cov
	5  144 cov
	6 132 cov
	RANDOM_RESIDUAL VALUES
	2.0

| **things to consider in the parameter template file**
| 1. Note that in the current version, the MLMs are tested with two mixed effects : additive polygenic effect and additive QTL effects. In case of asreml, these two mixed effects must be included in the parameter file and must be defined exactly with the same keywords (iDip and indi) as shown in the parameter file. 
| 2. In case of blupf90+, any additional option can be included after the last line (showing the Random residual values). 

