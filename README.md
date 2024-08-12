# Introduction

popgen-cldla is a fully automated workflow to identify QTLs using
combined linkage disequilibrium and linkage analysis (CLDLA). Unlike
genome-wide association study (GWAS), cLDLA is a mixed-linear-model
(MLM) based apporach to estimate variance component of each putative QTL
using the specialized tools such as ASREML and blupf90+.

The workflow consists of the following steps:

1. Phase the genotypes using beagle or shapeit (optional)
2. Filtering the SNPs (optional)
3. Estimation of genetic relationship across genome (GRM)
4. Estimation of genetic relationship at each putative QTL region
   (DRM)
5. Bending and inverting GRM and DRM
6. Variance component estimation
7. Plotting the results in the form of manhattan plot

**Note that the workflow currently only supports the singularity
container**

# Input files

1. List of phased or unphased vcf files in a csv format

```Bash
10,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr10.TxT.vcf.gz,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr10.TxT.vcf.gz.csi
11,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr11.TxT.vcf.gz,/data/testing/TailMLS04/OUT_VCF_BEAGLE4_ALL_OARAutoSom1_Chr11.TxT.vcf.gz.csi
.
.
.
```

2. Phenotype file

```Bash
iidx  iid  sex  age  bodywgt  witherhgt  taillgt
1  MLS169  1  11  10.6  47.0  32.7
2  MLS170  1  9  12.0  44.2  33.6
.
.
.
```

_Note that the header should be present in the phenotype file and
first, second and the last column in the phenotype file should be the
index of individuals (with the header: iidx), individual id (with the
header: iid) and the trait (Y), respectively._

3. Parameter file

Depending on the tool specified by the user, the workflow will
generate the parameter file.
**things to consider in the parameter template file**

1. If you have not used any of the two programs (asreml and blupf90+)
   implemented for variance component estimation before, please Refer to
   its respective documentation: [asreml
   documentation](https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-R-Reference-Manual-4.2.pdf)
   or [blupf90+
   documentation](http://nce.ads.uga.edu/html/projects/programs/docs/blupf90_all8.pdf).
   | 2. Note that in the current version, the MLMs are tested with two
   mixed effects : additive polygenic effect and additive QTL effects.

# Example commands to run the workflow

**To identify QTLs using CLDLA approach**

```Bash
nextflow run popgen-cldla/ --input chrom_vcf_idx.csv --maf 0.05 --pheno_file TailMLS04.template.phe -qs 10 --outdir testing_blupf90_TailMLS04 -resume -profile singularity --output_prefix TailMLS04 --tool blupf90
```

**To estimate heritability using the approach as implemented in GCTA**

```Bash
nextflow run popgen-cldla/ --input nextflow_testing/TailMLS04/chrom_vcf_idx.csv --maf 0.05 --pheno_file nextflow_testing/TailMLS04/TailMLS04.template.phe -qs 10 --outdir testing_h2_TailMLS04 -resume -profile singularity --output_prefix TailMLS04_h2 --estimate_h2
```

_Note that estimation of heritability using GCTA requires that the
regressors be separated into two files: quantitative variables and
qualitative variables. Therefore, in the phenotype file, any column with
the float values (identified using the presence of dot,\".\") are
automatically classified as quantitative and the column without float
values are classified as qualitative. Further, the workflow to estimate
h2, will produce the error if there is any column with mixtures of float
and integer values._

# Description of the parameters

```Bash
--input                       [string]  Path to comma-separated file containing information about the samples in the experiment.
--outdir                      [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                        infrastructure.
--email                       [string]  Email address for completion summary.
--maf                         [number]  minor allele frequency threshold; SNPs with MAF less than this threshold will filtered out
--window_size                 [number]  Window size to carry out cLDLA analysis [default: 40]
--output_prefix               [string]  output prefix should not contain the dot in it [default: cldla_run1]
--pheno_file                  [string]  path to the phenotypes file as recognized by echidna
--p_value                     [number]  p-value cutoff based on permutation test
--include_chrom               [string]  file containing ids of chromosome on which cLDLA will be carried out
--num_autosomes               [number]  total number of autosomes in the dataset
--estimate_h2                 [boolean] whether to estimate heritability using gcta tool
--lrt_threhold                [number]  lrt values above which the values are significant [default: 18]
--tool                        [string]  tool to estimate variance components [default: blupf90]
--par_file                    [string]  parameter file of the tool
--n_perm                      [integer] number of permutation test to be carried out to determine the significant threshold for CLDLA or H2 estimation
                            [default: 100]
--phase_genotypes             [boolean] whether or not to phase the genotypes
--phasing_panel               [string]  csv file containing information about path to the vcf files to be used for imputation
--phasing_map                 [string]  csv file containing information about path to the recombination map files
--phasing_tool                [string]  tool to be used for phasing: beagle5 or shapeit5 [default: beagle5]
```
