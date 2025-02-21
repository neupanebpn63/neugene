#' Run GWAS with GLM and Enhanced Visualization
#'
#' @description Performs GWAS using linear regression with improved visualization
#' @param y Phenotype data frame
#' @param X Genotype data frame
#' @param C Covariate data frame
#' @param snp_info SNP metadata data frame
#' @return List containing results and plots
#' @export
run_gwas_glm <- function(y, X, C, snp_info) {
  library(qqman)
  library(ggplot2)

  # Data preparation
  y_vec <- as.numeric(y[, 2])
  SNPs <- as.matrix(X[, -1])
  Covs <- as.matrix(C[, -1])

  # GWAS analysis
  p_values <- sapply(1:ncol(SNPs), function(i) {
    model <- lm(y_vec ~ SNPs[, i] + Covs)
    summary(model)$coefficients[2, 4]
  })

  # Prepare results
  results <- data.frame(SNP = colnames(SNPs), P = p_values) |>
    merge(snp_info, by = "SNP")

  # Enhanced Manhattan plot
  manhattan_plot <- function() {
    manhattan(results,
              chr = "Chromosome",
              bp = "Position",
              p = "P",
              snp = "SNP",
              col = c("#377eb8", "#4daf4a"),
              genomewideline = -log10(5e-8),
              suggestiveline = -log10(1e-5),
              main = "Manhattan Plot (GLM)",
              cex = 0.8,
              cex.axis = 0.9)
  }

  # Enhanced QQ plot
  qq_plot <- function() {
    qq(results$P,
       main = "QQ Plot (GLM)",
       col = "#984ea3",
       cex = 0.8)
    abline(0, 1, col = "#ff7f00", lwd = 2)
  }

  list(results = results,
       manhattan_plot = manhattan_plot,
       qq_plot = qq_plot)
}
