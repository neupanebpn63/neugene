% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gwas_glm_pca.R
\name{run_gwas_glm_pca}
\alias{run_gwas_glm_pca}
\title{Run GWAS with GLM + PCA Correction and Enhanced Visualization}
\usage{
run_gwas_glm_pca(y, X, C, snp_info, npc = 5)
}
\arguments{
\item{y}{Phenotype data frame}

\item{X}{Genotype data frame}

\item{C}{Covariate data frame}

\item{snp_info}{SNP metadata data frame}

\item{npc}{Number of principal components}
}
\value{
List containing results and plots
}
\description{
Performs PCA-adjusted GWAS with comprehensive visualization
}
