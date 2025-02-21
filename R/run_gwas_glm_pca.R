#' Run GWAS with GLM + PCA Correction and Enhanced Visualization
#'
#' @description Performs PCA-adjusted GWAS with comprehensive visualization
#' @inheritParams run_gwas_glm
#' @param npc Number of principal components
#' @return List containing results and plots
#' @export
run_gwas_glm_pca <- function(y, X, C, snp_info, npc = 5) {
  library(qqman)
  library(ggplot2)
  library(gridExtra)

  # Data preparation
  y_vec <- as.numeric(y[, 2])
  X_snps <- as.matrix(X[, -1])
  Covs <- as.matrix(C[, -1])

  # PCA analysis
  pca_out <- prcomp(X_snps, center = TRUE, scale. = TRUE)
  PCs <- pca_out$x[, 1:npc]

  # Collinearity check
  design_matrix <- Covs
  current_rank <- qr(design_matrix)$rank
  valid_pcs <- numeric(0)

  for(i in 1:npc) {
    temp_design <- cbind(design_matrix, PCs[, i])
    if(qr(temp_design)$rank > current_rank) {
      valid_pcs <- c(valid_pcs, i)
      design_matrix <- temp_design
      current_rank <- qr(temp_design)$rank
    }
  }

  # GWAS with PCA
  C_adj <- if(length(valid_pcs) > 0) cbind(Covs, PCs[, valid_pcs]) else Covs
  p_values <- sapply(1:ncol(X_snps), function(i) {
    model <- lm(y_vec ~ X_snps[, i] + C_adj)
    summary(model)$coefficients[2, 4]
  })

  # Prepare results
  results_pca <- data.frame(SNP = colnames(X_snps), P = p_values) |>
    merge(snp_info, by = "SNP")

  # Visualization functions
  manhattan_plot <- function() {
    manhattan(results_pca,
              chr = "Chromosome",
              bp = "Position",
              p = "P",
              snp = "SNP",
              col = c("#e41a1c", "#377eb8"),
              genomewideline = -log10(5e-8),
              suggestiveline = -log10(1e-5),
              main = "Manhattan Plot (GLM+PCA)",
              cex = 0.8,
              cex.axis = 0.9)
  }

  qq_plot <- function() {
    qq(results_pca$P,
       main = "QQ Plot (GLM+PCA)",
       col = "#984ea3",
       cex = 0.8)
    abline(0, 1, col = "#ff7f00", lwd = 2)
  }

  # PCA visualizations
  variance <- pca_out$sdev^2 / sum(pca_out$sdev^2)

  scree_plot <- ggplot(data.frame(PC = 1:npc, Variance = variance[1:npc]),
                       aes(x = PC, y = Variance)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_line(color = "darkred", size = 1) +
    geom_point(size = 3, color = "darkred") +
    labs(title = "PCA Scree Plot", x = "Principal Component", y = "Variance Explained") +
    theme_minimal()

  pc_scatter <- ggplot(data.frame(PC1 = pca_out$x[,1], PC2 = pca_out$x[,2]),
                       aes(x = PC1, y = PC2)) +
    geom_point(color = "#4daf4a", alpha = 0.6) +
    labs(title = "PC1 vs PC2", x = "Principal Component 1", y = "Principal Component 2") +
    theme_minimal()

  list(results = results_pca,
       manhattan_plot = manhattan_plot,
       qq_plot = qq_plot,
       pca_plots = list(scree_plot, pc_scatter),
       PCs_used = valid_pcs)
}
