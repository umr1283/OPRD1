Effectif_mutation_YCC <- function(Phenotypes, G_mat, y_name) {
  if (ncol(G_mat) <= 1) {
    Effectifs_CC <- data.frame(matrix(
      data = NA, 
      nrow = 1, 
      ncol = 9, 
      dimnames = list(
        NULL, 
        c(
          "Variant", 
          map_chr(.x = transpose(expand.grid(
            stringsAsFactors = FALSE,
            c("genotype_Missing", "genotype_0", "genotype_1", "genotype_2"),
            c("CTRL", "CASE")
          )), .f = ~paste(., collapse = "_"))
        )
      )
    ))[-1, ]
    return(Effectifs_CC)
  }
  
  counts_table <- map(.x = c(0, 1), .f = function(.x) {
    samples <- Phenotypes[Phenotypes[, y_name] %in% .x, "IID"]
    
    Mat_eff <- G_mat %>% 
      filter(IID %in% samples) %>%
      select(-IID) %>% 
      Zoom_Genotype() %>%
      select(Variant, Missing_genotype, G_0, G_1, G_2)
    
    Mat_percent <- Mat_eff %>% 
      transpose() %>% 
      map_df(.f = function(.l) {
        out <- unlist(.l[-grep("Variant", names(.l))])
        as_tibble(matrix(
          data = percent(out / sum(out, na.rm = TRUE)), 
          byrow = TRUE, 
          nrow = 1, dimnames = list(NULL, names(out))
        )) %>% 
          mutate(Variant = .l$Variant)
      }) %>% 
      dplyr::rename(
       Missing_genotype_percent = Missing_genotype,
       genotype_0_percent = G_0,
       genotype_1_percent = G_1,
       genotype_2_percent = G_2
      )
    
    percent_Eff <- inner_join(x = Mat_percent, y = Mat_eff, by = "Variant") %>% 
      mutate(
        genotype_Missing = paste0(Missing_genotype_percent, " (", Missing_genotype, ")"),
        genotype_0 = paste0(genotype_0_percent, " (", G_0, ")"),
        genotype_1 = paste0(genotype_1_percent, " (", G_1, ")"),
        genotype_2 = paste0(genotype_2_percent, " (", G_2, ")")
      ) %>% 
      select(Variant, genotype_Missing, genotype_0, genotype_1, genotype_2)
    
    new_colnames <- switch(EXPR = as.character(.x),
      "1" = colnames(percent_Eff)[-1] <- paste0(colnames(percent_Eff)[-1], "_CASE"),
      "0" = colnames(percent_Eff)[-1] <- paste0(colnames(percent_Eff)[-1], "_CTRL")
    )
    colnames(percent_Eff)[-1] <- new_colnames
    return(percent_Eff)
  })
  
  Effectifs_CC <- inner_join(x = counts_table[[1]], y = counts_table[[2]], by = "Variant")
  return(Effectifs_CC)
}
