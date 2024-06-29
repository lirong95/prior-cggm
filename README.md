# prior-cggm
Incorporating prior information in gene expression network-based cancer heterogeneity analysis

1. functions.R: contains all the functions in the proposed algorithm. Please source it directly for use.
2. metabric_brca_pampath.RData: cleaned breast cancer dataset from METABRIC.
   covariate (x): copy number variation matrix
   data (y): mRNA gene expression matrix
   prior_omega: prior infromation matrix extracted by easyPubMed
3. metabric_brca_example_code.R: the code in the real data analysis.
