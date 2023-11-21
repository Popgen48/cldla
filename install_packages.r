if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#impute_url<-"https://bioconductor.org/packages/release/bioc/src/contrib/impute_1.76.0.tar.gz"
#install.packages(impute_url, repos=NULL, type="source")
#
#snpready_url<-"https://cran.r-project.org/src/contrib/snpReady_0.9.6.tar.gz"
#install.packages(snpready_url, repos=NULL, type="source")
#
#aghmatrix_url<-"https://cran.r-project.org/src/contrib/AGHmatrix_2.1.4.tar.gz"
#install.packages(aghmatrix_url, repos=NULL, type="source")

install.packages("AGHmatrix")
BiocManager::install("impute")

install.packages("snpReady")
