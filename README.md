# kellets_whelk_rnaseq

## Project description
This is a project analyzing population structure and genetic adaptation of Kellet's whelks (Kelletia kelletii)  --- a range expanding marine species with high gene flow potential. We performed a series of crosses on Kellet's whelks collected from its historical and recently colonized range, and conducted RNA-Seq on 70 offspring that we reared in a common garden environment. We also used transcriptome derived SNPs to examine population structure and genetic adaptation. 

Corresponding article [place holder for citation]

## Scripts
Code for project is in [scripts](https://github.com/ChristieLab/kellets_whelk_rnaseq/tree/main/scripts). Analyses are split up into directories: 
1. _**Study system.**_ R code associated with figure 1.
2. _**Sequence processing.**_ Process RNA transcripts and conduct DEG analysis
3. _**SNP calling.**_ Contains for SNP variant calling using ANGSD 
4. _**Population structure and population assignment.**_ R codes analyszing transcriptome-derived SNPs using snpR (Figure 2)
5. _**Candidate loci driving differntiation.**_  R codes analyszing transcriptome-derived SNPs using snpR, population assignment using rubias and ranger (Figure 3). R code analyzing the TPI gene (Figure 4)
6. _**Gene Ontology and WGCNA.**_ R codes associated with analyses of expression data (Supplemental)  
  

## Data
All raw data and sample metadata generated for this project are stored in the NCBI Short Read Archive (SRA) under project PRJNA1000198 (https://www.ncbi.nlm.nih.gov/sra/PRJNA1000198). VCF files for the all-SNPs and DEG-SNPs datasets are deposited in Dryad (https://doi.org/10.5061/dryad.qbzkh18s3). (https://github.com/ChristieLab/kellets_whelk_rnaseq/assets/86313796/02b795ee-c350-43cf-bb37-8bee7a0f3e9b)
