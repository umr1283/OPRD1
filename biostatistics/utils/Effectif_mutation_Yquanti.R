Effectif_mutation_Yquanti <- function(G_mat) {
  if (ncol(G_mat) == 0) {
    Pourc_Eff <- data.frame(matrix(rep(NA, 4), byrow = TRUE, nrow = 1))[-1, ]
    names(Pourc_Eff) <- c("Variant", "genotype_0", "genotype_1", "genotype_2")
  } else {
    detail_analysis_real <- Zoom_Genotype(G_mat) %>% select("Variant", "G_0", "G_1", "G_2")

    Mat_pourc <- round(t(apply(data.matrix(detail_analysis_real[, -1]), 1, prop.table)) * 100, 2)
    Mat_pourc <- cbind.data.frame("Variant" = detail_analysis_real$Variant, Mat_pourc)
    names(Mat_pourc) <- c("Variant", "G_0_pourc", "G_1_pourc", "G_2_pourc")

    Pourc_Eff <- merge(Mat_pourc, detail_analysis_real, by = "Variant")

    Pourc_Eff$genotype_0 <- paste0(Pourc_Eff$G_0_pourc, "% (", Pourc_Eff$G_0, ")")
    Pourc_Eff$genotype_1 <- paste0(Pourc_Eff$G_1_pourc, "% (", Pourc_Eff$G_1, ")")
    Pourc_Eff$genotype_2 <- paste0(Pourc_Eff$G_2_pourc, "% (", Pourc_Eff$G_2, ")")

    Pourc_Eff <- Pourc_Eff[, c("Variant", "genotype_0", "genotype_1", "genotype_2")]
  }
  return(Pourc_Eff)
}
