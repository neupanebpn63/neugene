
# neugene: Genome-Wide Association Analysis Using GLM and PCA Correction

**neugene** is an R package for conducting genome-wide association studies (GWAS) using general linear models (GLM), with options for incorporating principal component analysis (PCA) to account for population structure.

The package provides SNP-level association testing, optional PCA-based covariate adjustment, and visualization of GWAS results through Manhattan plots, quantile–quantile (QQ) plots, and PCA plots. It is designed to provide a straightforward and accessible workflow for researchers working with genotype and phenotype datasets.

The package is developed and maintained by **Bipin Neupane** at Washington State University.

## Features

- **GLM-based GWAS:** Performs single-marker association testing using linear regression.
- **PCA-corrected GWAS:** Incorporates principal components derived from genotype data as additional covariates.
- **Collinearity assessment:** Identifies principal components that contribute additional independent information relative to the supplied covariates.
- **Manhattan plots:** Visualizes genome-wide SNP association signals across chromosomes.
- **QQ plots:** Compares observed and expected association P-value distributions.
- **PCA visualization:** Generates scree plots and PC1–PC2 scatter plots.
- **Exportable results:** Returns SNP association results, plotting functions, and information about retained principal components.

## Installation

The development version of `neugene` can be installed directly from GitHub using the `remotes` package.

```r
# Install remotes if needed
install.packages("remotes")

# Install neugene from GitHub
remotes::install_github("neupanebpn63/neugene")

# Load the package
library(neugene)
```

The package uses `ggplot2`, `qqman`, and `rlang` for visualization and supporting functionality. Required dependencies are declared in the package's DESCRIPTION file.

## Input Data

Both GWAS functions require four input datasets:

| Argument | Description |
|----------|-------------|
| `y` | Phenotype data frame with sample identifiers in the first column and numeric phenotype values in the second column. |
| `X` | Genotype data frame with sample identifiers in the first column and SNP genotype scores in subsequent columns. |
| `C` | Covariate data frame with sample identifiers in the first column and numeric covariates in subsequent columns. |
| `snp_info` | SNP metadata containing `SNP`, `Chromosome`, and `Position` columns. |

**Important:** The phenotype, genotype, and covariate data must contain the same individuals in the same order. The current functions assume that sample alignment has been performed before analysis.

SNP genotype columns should be numeric, and SNP identifiers should match those in `snp_info`. Missing values, monomorphic markers, and incomplete metadata should be addressed during data preparation.

## Usage

### 1. GWAS Using General Linear Model (GLM)

The `run_gwas_glm()` function performs single-marker association testing using linear regression, incorporating the supplied covariates.

```r
library(neugene)

# y, X, C, and snp_info are preprocessed input data frames

gwas_glm <- run_gwas_glm(
  y = y,
  X = X,
  C = C,
  snp_info = snp_info
)

# View association results
head(gwas_glm$results)

# Generate Manhattan plot
gwas_glm$manhattan_plot()

# Generate QQ plot
gwas_glm$qq_plot()
```

The returned association table includes SNP identifiers, P-values, and corresponding marker metadata.

### 2. GWAS with PCA Correction

The `run_gwas_glm_pca()` function extends the GLM approach by incorporating principal components calculated from the genotype matrix.

The number of principal components can be specified using the `npc` argument (default = 5).

```r
gwas_pca <- run_gwas_glm_pca(
  y = y,
  X = X,
  C = C,
  snp_info = snp_info,
  npc = 5
)

# View association results
head(gwas_pca$results)

# Generate Manhattan plot
gwas_pca$manhattan_plot()

# Generate QQ plot
gwas_pca$qq_plot()

# PCA scree plot
print(gwas_pca$pca_plots[[1]])

# PC1 vs PC2 scatter plot
print(gwas_pca$pca_plots[[2]])

# View retained principal components
gwas_pca$PCs_used
```

The PCA-corrected model evaluates whether each principal component adds independent information to the existing covariate matrix. Components that do not increase the matrix rank are excluded from the final covariate set.

### 3. Exporting GWAS Results

Association results can be saved as CSV files for further analysis.

```r
# Export GLM results
write.csv(
  gwas_glm$results,
  "GLM_GWAS_Results.csv",
  row.names = FALSE
)

# Export PCA-corrected GLM results
write.csv(
  gwas_pca$results,
  "GLM_PCA_GWAS_Results.csv",
  row.names = FALSE
)
```

## Example Datasets and Tutorials

The repository includes example genotype, phenotype, covariate, and SNP metadata files for exploring the package's functionality.

**Example input datasets:**

- `Gapit_genotype.txt` — genotype data
- `phenotype.txt` — phenotype data
- `covariates.txt` — covariate data
- `snp_info.txt` — SNP chromosome and position information

These datasets are organized under `inst/extdata/` and can be located after installation using:

```r
system.file(
  "extdata",
  "phenotype.txt",
  package = "neugene"
)
```

**Additional resources:**

- [HTML Tutorial](examples/neugene_tutorial.html)
- [PDF Tutorial](examples/neugene_tutorial.pdf)
- [GLM GWAS Results](examples/GLM_GWAS_Results.xlsx)
- [PCA-Corrected GWAS Results](examples/GLM_PCA_GWAS_Results.xlsx)
- [GWAS Results and Visualizations](examples/GWAS_Results.pdf)

The tutorials and example outputs provide supplementary guidance for understanding the GWAS workflow and interpreting association results.

## Methodological Considerations

The package implements single-marker GLM association testing, with optional PCA adjustment for population structure.

The following considerations should be kept in mind:

- Sample identifiers and row order must be consistent across phenotype, genotype, and covariate datasets.
- Genotype matrices used for PCA should contain suitable numeric marker scores without missing or zero-variance columns.
- The number of requested PCs should be appropriate for the dimensions of the genotype dataset.
- The current implementation uses fixed Manhattan plot thresholds of 5 × 10⁻⁸ (genome-wide) and 1 × 10⁻⁵ (suggestive). These thresholds should not be assumed to be appropriate for every crop or GWAS design.
- The package does not implement a mixed linear model or genomic relationship matrix correction. PCA adjustment alone may not account for all confounding due to relatedness or population structure.
- The package assumes SNP regression coefficients are estimable. Additional preprocessing may be required for problematic or monomorphic markers.

Users should evaluate model assumptions, marker quality, population structure, and multiple-testing thresholds according to their research objectives.

## Package Development and Validation

`neugene` is developed as an R package with documented functions, declared dependencies, and an MIT license.

The package has been evaluated using `devtools::check()` on Windows with R 4.4.3, obtaining zero errors and zero warnings. This confirms package-level structural and dependency checks, but does not independently establish the statistical validity of every analysis or input scenario.

The package is under active development, and future updates may introduce additional functionality, input validation, and improvements in computational efficiency.

## Citation

If you use `neugene` in your research, please cite the software repository and the versioned Zenodo archive.

**Software:** Neupane, B. (2026). *neugene: Genome-Wide Association Analysis Using GLM and PCA Correction*. R package.

**GitHub:** https://github.com/neupanebpn63/neugene

**Zenodo DOI:** To be added after archival.

A permanent, version-specific DOI will be provided through Zenodo to preserve the software version associated with research applications.

## License

This package is distributed under the [MIT License](LICENSE.md).

## Author and Contact

**Bipin Neupane**  
Washington State University  
Pullman, Washington, USA

Email: neupanebpn63@gmail.com  
GitHub: [neupanebpn63](https://github.com/neupanebpn63)
Website: [neupanebpn](https://neupanebpn63.wixsite.com/bpn1)

For questions, bug reports, or suggestions, please open an issue in the [GitHub repository](https://github.com/neupanebpn63/neugene/issues).
