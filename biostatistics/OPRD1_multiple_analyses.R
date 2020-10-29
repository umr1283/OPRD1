#!/usr/bin/env Rscript

################################
#### MAKE ANALYSES OF OPRD1 ####
################################

## ex : 
## nohup Rscript OPRD1_multiple_analyses.R GPCR GOF CC_OBESITE TRUE & 
## nohup Rscript OPRD1_multiple_analyses.R GPCR LOF CC_OBESITE TRUE & 
## nohup Rscript OPRD1_multiple_analyses.R GPCR ALL CC_OBESITE TRUE & 
## nohup Rscript OPRD1_multiple_analyses.R GPCR LOFgene+LOFfunc CC_OBESITE TRUE & 


options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

project_directory <- "/disks/PROJECT/OPRD1"

#### Load packages ---------------------------------------------------------------------------------
suppressPackageStartupMessages(
  invisible(
    lapply(
      c(
        "datasets", "utils", "grDevices", "graphics", "stats", "methods",
        "parallel", "scales", "readxl", "writexl", "tidyverse", "mongolite"
      ),
      library,
      character.only = TRUE
    )
  )
)

RawDir <- "/sequencing_data/"
PhenotypeDir <- "/phenotypes_data/"
Phenotype_file <- paste0(PhenotypeDir, "Phenotyes_latest.xlsx")
SourceDir <- paste0(project_directory, "/Scripts/utils")
InputDir <- paste0(project_directory, "/Data/00_MAKE_DATASET/")
OutDir <- paste0(project_directory, "/Data/OPRD1_multiple_analyses/")

## to adapt multiple analyses ##
Pheno_OPRD1 <- "/disks/PROJECT/OPRD1/Docs/OPRD1-pheno.xlsx" 
Analyses_plan <- read_excel(
  path = "/disks/PROJECT/OPRD1/Docs/Analyses_plan_OPRD1.xlsx", 
  col_names = TRUE,
  na = "NA",
) %>%
  as.data.frame()

#### Functions -------------------------------------------------------------------------------------
source(paste0(SourceDir, "/MiST.R"))
source(paste0(SourceDir, "/Parameters-class.R"))
source(paste0(SourceDir, "/QC-class.R"))
source(paste0(SourceDir, "/Summary_clinique.R"))
source(paste0(SourceDir, "/Zoom_genotype.R"))
source(paste0(SourceDir, "/Effectif_mutation_YCC.R"))
source(paste0(SourceDir, "/Effectif_mutation_Yquanti.R"))
source(paste0(SourceDir, "/DO_QC.R"))

f <- function(objectParams) {
  slot_names <- slotNames(objectParams)
  out <- as.list(rep(NA, length(slot_names))) %>%
    `names<-`(slot_names)
  for (islot in slot_names) {
    out[[islot]] <- objectParams[islot]
  }
  return(as.data.frame(out))
}
f2 <- function(objectQC) {
  out <- data.frame(matrix(objectQC["Excluded"], byrow = TRUE, ncol = 1)) %>%
    `names<-`(paste0("QC_", objectQC["Step"], "_", objectQC["Type"]))
  return(out)
}
f3 <- function(listf2) {
  nrowQCtable <- max(unlist(lapply(listf2, nrow)))
  out <- data.frame(matrix(nrow = nrowQCtable, ncol = length(listf2)))
  names(out) <- unlist(lapply(listf2, names))
  for (el in listf2) {
    if (nrow(el) < nrowQCtable) {
      out[, names(el)] <- c(el[, names(el)], rep(NA, nrowQCtable - nrow(el)))
    } else {
      out[, names(el)] <- el
    }
  }
  return(out)
}
f4 <- function(objectQC) {
  listout <- list(
    "Why" = objectQC["Why"],
    "Type" = objectQC["Type"],
    "N_origin" = objectQC["N_origin"],
    "Excluded" = objectQC["Excluded"],
    "N_excluded" = objectQC["N_excluded"],
    "Comment" = objectQC["Comment"]
  )
  return(listout)
}

#### DEF -------------------------------------------------------------------------------------------
LARGEdefSNP <- c(
  "initiator_codon_variant", "synonymous_variant", "stop_retained_variant", "stop_lost", "stop_gained",
  "missense_variant", "incomplete_terminal_codon_variant", "coding_sequence_variant", "splice_acceptor_variant",
  "splice_donor_variant", "splice_region_variant", "regulatory_region_variant"
)
STRICTdefSNP <- c(
  "initiator_codon_variant", "stop_retained_variant", "stop_lost", "stop_gained", "missense_variant",
  "splice_donor_variant", "splice_acceptor_variant"
)
LARGEdefINDEL <- c("all")
STRICTdefINDEL <- c("frameshift_variant", "inframe_insertion", "splice_acceptor_variant", "splice_donor_variant")

#### SET PARAMETERS --------------------------------------------------------------------------------
switch(
  EXPR = length(args),
  "1" = {
    args[2] <- "LARGE"
    args[3] <- "CC_T2D61"
    args[4] <- "TRUE"
    message(paste("[MONO-GENE]", "The group", args[1], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Analysis LARGE is selected, by default!"))
    message(paste("[MONO-GENE]", "Pheno CC_T2D61 is selected, by default!"))
    message(paste("[MONO-GENE]", "Pheno is set BINARY = TRUE"))
  },
  "2" = {
    message(paste("[MONO-GENE]", "The group", args[1], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Analysis", args[2], "was found in 'args'!"))
    args[3] <- "CC_T2D61"
    message(paste("[MONO-GENE]", "Pheno CC_T2D61 is selected, by default!"))
    args[4] <- "TRUE"
    message(paste("[MONO-GENE]", "Pheno is set BINARY = TRUE"))
  },
  "3" = {
    message(paste("[MONO-GENE]", "The group", args[1], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Analysis", args[2], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Pheno", args[3], "was found in 'args'!"))
    args[4] <- "TRUE"
    message(paste("[MONO-GENE]", "Pheno is set BINARY = TRUE"))
  },
  "4" = {
    message(paste("[MONO-GENE]", "The group", args[1], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Analysis", args[2], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Pheno", args[3], "was found in 'args'!"))
    message(paste("[MONO-GENE]", "Pheno is set BINARY =", args[4]))
  },
  stop(
    paste(
      "[MONO-GENE]",
      "4 arguments must be supplied GroupeGene, LARGE_STRICT, the pheno studied, is it BINARY or NOT (TRUE/FALSE).",
    ),
    call. = FALSE
  )
)

GroupeGene <- as.character(args[1])
analyse_LARGE_STRICT <- as.character(args[2])
y_name <- as.character(args[3])
binary <- as.logical(toupper(as.character(args[4])))

## SET PARAM 
SNP_INDEL <- "BOTH"
UseCustomFile <- FALSE
CustomFile <- paste0(project_directory, "/Docs/CUSTOM_files/_nocustom_toDO.xlsx")
PopFilter <- "ALL"

Gene_to_run <- "OPRD1"
ENST_to_run <- "ENST00000234961"


#### Selection of genes and connexion  -------------------------------------------------------------
liste_TR <- read.table(
  file = paste0(RawDir, "Docs/Latest_info_gene.csv"),
  header = TRUE,
  sep = ";"
) %>%
  filter(get(GroupeGene) == 1) %>%
  select(gene_name, Transcript) %>%
  filter(gene_name %in% Gene_to_run)

message(paste("[MONO-GENE]", "<===>"))
message(paste("[DATE-TIME]",Sys.time()))
message(paste("[MONO-GENE]", "Group selected:", GroupeGene))
message(paste("[MONO-GENE]", "Analysis:", analyse_LARGE_STRICT))
message(paste("[MONO-GENE]", "Pheno:", y_name))
message(paste("[MONO-GENE]", "Binary:", binary))

## read Pheno
Phenotypes <- read_excel(
  path = Phenotype_file,
  col_names = TRUE,
  na = "NA",
  guess_max = 9365
) %>%
  as.data.frame()
rownames(Phenotypes) <- paste(Phenotypes$RUN, Phenotypes$ID, sep = "_")
  
#### Load specifiq phenotype for OPRD1 abalyses -------------------------------------------------
Phenotypes_OPRD1 <- read_excel(
  path = Pheno_OPRD1, 
  col_names = TRUE,
  na = "NA",
  guess_max = 9365
) %>%
  as.data.frame() %>% 
  mutate(ADN = as.character(ADN))

#### Update Phenotypes /!\ 
Phenotypes <- Phenotypes %>% 
  filter(ADN %in% Phenotypes_OPRD1$ADN) %>% 
  select(-names(Phenotypes_OPRD1)[names(Phenotypes_OPRD1)%in% names(Phenotypes)][-1]) %>% 
  left_join(x = ., y = Phenotypes_OPRD1, by = "ADN") %>% 
  mutate(
    age = AGE_D0, 
    sex = SEX, ## must be a numeric in the MiST matrix 
    CC_OWT = `CC-OWT`, 
    CC_OBESITE_avecINFODIAB = CC_OBESITE, 
    FG = FG_D0, 
    FI = FI_D0, 
    LN_FI = log(x = FI),
    HOMA2B = HOMA2B_D0, 
    LN_HOMA2B = log(x = HOMA2B_D0), 
    HOMA2IR = HOMA2IR_D0,
    LN_HOMA2IR = log(x = HOMA2IR_D0), 
    HDL = `HDL-sansTTT`, 
    TG = `TG-sansTTT`, 
    LN_TG = log(x = `TG-sansTTT`), 
    PAS = `PAS-sansTTT`, 
    PAD = `PAD-sansTTT`, 
    CC_MS_sansBMI = CC_MS, 
    LN_BMI = log(BMI)
  ) 


## make plot => histogram or density plot (with transformation) unique(Analyses_plan$Y)
if (FALSE) {
  tmp <- Phenotypes %>%
    mutate_at(
      .vars = vars(unique(Analyses_plan$Y[Analyses_plan$Type %in% "binaire"])),
      .funs = function(.x) as.factor(as.character(.x))
    )

  y_interest <- unique(Analyses_plan$Y)
  y_interest <- y_interest[-grep("_avec|_sans", y_interest)]
  y_trans <- y_interest[grep("LN_", y_interest)]
  y_interest <- y_interest[!y_interest %in% y_trans]
  y_interest <- c(
    sort(y_interest), 
    c(matrix(data = c(y_trans, gsub("LN_", "", y_trans)), nrow = 2, byrow = TRUE))
  )
  numeric_var <- names(dplyr::select_if(tmp, is.numeric))
  for (ivar in numeric_var[numeric_var %in% y_interest]) {
    cat(paste("##", ivar, " \n"))
    (ggplot(data = tmp, mapping = aes(x = get(ivar))) + 
      geom_density(aes(y = ..count..)) +
      labs(x = ivar)
    ) %>% print()
    cat("\n\n")
  } 
  fac_var <- names(dplyr::select_if(tmp, is.factor))
  for (ivar in fac_var[fac_var %in% y_interest]) {
    cat(paste("##", ivar, " \n"))
    (ggplot(data = tmp, mapping = aes(x = get(ivar))) + 
      geom_histogram(stat="count") +
      labs(x = ivar)
    ) %>% print()
    cat("\n\n")
  } 
}

#### BEGIN OF ANALYSIS ---------------------------------------------------------------------------
options(warn = 2) ## Turn warning into error => can be trapped

{

  #### ALL parameters ----------------------------------------------------------------------------
  Params <- new.OPRD1_Parameters()
  Params["FolderOut"] <- OutDir
  Params["GroupeGene"] <- GroupeGene
  Params["Gene_Symbol"] <- Gene_to_run
  Params["ENST"] <- ENST_to_run
  Params["analyse_LARGE_STRICT"] <- analyse_LARGE_STRICT
  Params["threshold"] <- 0.01
  Params["ExludeSynonymous"] <- FALSE
  Params["cluster"] <- "none"
  Params["y_name"] <- y_name
  Params["binary"] <- binary
  Params["SNP_INDEL"] <- SNP_INDEL
  # Params
  message(paste(
    "[MONO-GENE]", ">>> Go on", Params["Gene_Symbol"], Params["ENST"],
    "Version", Params["analyse_LARGE_STRICT"]
  ))
  
  dir.create(
    path = paste0(Params["FolderOut"], Params["Gene_Symbol"]),
    recursive = TRUE,
    showWarnings = FALSE,
    mode = "0777"
  )
  message(paste("[MONO-GENE]", "<Gene & Transcript ever seen, can load Phenotypes and RData..."))


  #### RData ready to analyses
  dataset_file <- paste0(
    InputDir,
    Params["Gene_Symbol"], "/", Params["Gene_Symbol"], "_", Params["ENST"], ".RData"
  )
  load(dataset_file)
  #### Update Genotypes /!\ with genotype listed in Pheno_OPRD1 => confirmed by SANGER -------------

  all_variants <- unique(c(
    unique(Phenotypes$LOF)[-1], 
    unique(Phenotypes$GOF)[-1], 
    unique(Phenotypes$ALL)[-1], 
    unique(Phenotypes$`LOFgene+LOFfunc`)[-1]
  ))
  
  all_genotype_data <- all_genotype_data[Phenotypes$IID, unique(Phenotypes[,Params["analyse_LARGE_STRICT"]])[-1]]
  all_annotation_data <- all_annotation_data %>% 
    filter(Variant %in% names(all_genotype_data))
  Detail_Geno <- Detail_Geno %>% 
    filter(Variant %in% names(all_genotype_data))
  all_quality_data <- all_quality_data %>% 
    filter(chr_pos %in% all_annotation_data$chr_pos)
  
  
  #### Select covariates -----------------------------------------------------------------------
  covar_names <- unlist(strsplit(
    x = unique(Analyses_plan[Analyses_plan$Y %in% Params["y_name"], "Covars"]), split = " + ", fixed = TRUE
  )) %>% gsub(" ", "", .)
  
  message("## Analysis can be done : Go !  ##")
 
  #### RUN ANALYSIS ------------------------------------------------------------------------------
  message(paste("[MONO-GENE]", "Go to run analysis...")) 
  
  #### KEEP ONLY SAMPLES INVOLVED IN THE STUDY ---------------------------------------------------
  message(paste("[MONO-GENE]", nrow(Phenotypes), "samples at the begining."))
  QC_PHENO <- Phenotypes %>%
    select("IID", Params["y_name"]) %>%
    `colnames<-`(c("IID", "y_name")) %>%
    filter(!is.na(y_name))
  
  ## Remove sample
  Phenotypes <- Phenotypes %>%
    filter(IID %in% QC_PHENO$IID)
  
  all_genotype_data <- all_genotype_data[Phenotypes$IID, , drop = FALSE]
  all_genotype_data$IID <- rownames(all_genotype_data)
  all_genotype_data <- all_genotype_data %>%
    filter(IID %in% QC_PHENO$IID)
  
  ## To study FG, only keep normoglycemic samples
  ## rule : FG < 7 and CC_Diab != 1
  if (Params["y_name"] %in% c("FG", "LN_FG")) {
    message(paste("[MONO-GENE]", "The study of FG involved to remove cases diabetique and child."))
    message(paste("[MONO-GENE]", "Keep only 'adult with FG<7 & CC_T2D56 %in% c(0,NA)'"))
    Phenotypes <- Phenotypes %>%
      filter(FG<7 & CC_T2D56 %in% c(0,NA) & selectInd %in% "adult")
    all_genotype_data <- all_genotype_data %>%
      filter(IID %in% Phenotypes$IID)
  }
  
  ## if option PopFilter == "EUR" ## 
  if (PopFilter %in% "EUR") {
    message("[MONO-GENE]", "Only EUR samples will be kept...")
    Phenotypes <- Phenotypes %>%
      filter(Pop_pred %in% "EUR")
    all_genotype_data <- all_genotype_data %>%
      filter(IID %in% Phenotypes$IID)
  }
  
  message(paste("[MONO-GENE]", nrow(Phenotypes), "samples involved in this study."))
  
  message(paste("[MONO-GENE]", ncol(all_genotype_data) - 1, "variants ready to be analysed."))
  
  #### QC steps ------------------------------------------------------------------------------
  QC <- DO_QC() ## all object (all_genotype, Phenotypes, ...) are updated. 
  QCdone <- TRUE
    
  if (ncol(all_genotype_data) == 1) {
    all_genotype_data <- data.frame(matrix(nrow = 0, ncol = 0))
    Phenotypes <- data.frame(matrix(nrow = 0, ncol = 0))
  }
  
  #### clinicalData --------------------------------------------------------------------------
  if (ncol(Phenotypes) != 0 & nrow(Phenotypes) != 0) {
    if (Params["binary"]) {
      clinicalData <- Summary_clinique(data = Phenotypes, x = Params["y_name"], na_symbol = "NA")
    } else {
      clinicalData <- summary(Phenotypes[, Params["y_name"]]) %>%
        t() %>%
        as.data.frame() %>%
        select(-Var1) %>%
        `colnames<-`(c(Params["y_name"], "Values")) %>%
        mutate(
          Values = round(x = Values, digit = 4)
        )
      clinicalData[, 1] <- as.character(clinicalData[, 1])
      clinicalData <- rbind.data.frame(clinicalData, c("N", nrow(Phenotypes)))
    }
    clinicalData <- list(clinicalData)
    message(paste("[MONO-GENE]", "clinicalData Done"))
  } else {
    clinicalData <- NULL
    message(paste("[MONO-GENE]", "/!\\ No phenotypes"))
  }
  
  #### Analysis using MiST -------------------------------------------------------------------
  message(paste("[MONO-GENE]", "Analysis using MiST"))
  variantRARE <- Detail_Geno$Variant[Detail_Geno$MAF < Params["threshold"]]
  G_rare <- all_genotype_data[, names(all_genotype_data) %in% variantRARE, drop = FALSE]
  message(paste("[MONO-GENE]", "G_rare Done"))
  if (ncol(G_rare) > 1) { ## rare analysis possible
    message(paste("[MONO-GENE]", "After all QCs,", ncol(G_rare), "rare variants will be analysed."))
    
    # param cluster = "none" ## or "TYPE"
    if (Params["cluster"] == "TYPE") {
      info_filtred <- Detail_Geno[Detail_Geno$Variant %in% names(G_rare), ]
      if (length(unique(info_filtred$Consequence)) > 1) {
        Z <- model.matrix(~ Consequence - 1, data = info_filtred) ## Z mat indicatrice
      } else { ## 1 cluster because 1 Type of consequence
        Z <- matrix(rep(1, ncol(G_rare)), ncol = 1)
      }
    message(paste("[MONO-GENE]", "WARNING cluster=TYPE selected, out_mist must be adapted..."))
    
    } else { ## cluster = "none" => 1 cluster
      Z <- matrix(rep(1, ncol(G_rare)), ncol = 1)
    }
    y <- Phenotypes[, Params["y_name"]]
    X <- as.matrix(Phenotypes[, covar_names, drop = FALSE])
    G <- as.matrix(G_rare)
    Z <- as.matrix(Z)
    GZ <- G %*% Z
      
    if (Params["binary"]) {
      
      ## mist ## 
      has_issues <- try(tools::assertCondition(out_mist <- mist(
        y = y, 
        X = X, 
        G = G, 
        Z = Z, 
        model = "binary"
      )), silent = TRUE)
      has_issues <- has_issues[
        sapply(has_issues, function(el) {length(intersect(class(el), c("warning", "error")))!=0})
      ]
      if (length(has_issues)==0) {
        out <- tibble::tibble(
          trait = Params["y_name"], 
          covariates = paste(covar_names, collapse = ";"),
          sample_size = nrow(Phenotypes),
          error = "None",
          mist_statistics = list(mist_print(out_mist)$statistic),
          mist_estimate = list(mist_print(out_mist)$estimate)
        ) 
      } else {
        out <- tibble::tibble(
          trait = Params["y_name"], 
          covariates = paste(covar_names, collapse = ";"),
          sample_size = nrow(Phenotypes),
          error = paste(
            unique(sapply(
              X = has_issues, 
              FUN = function(el) {paste0("[", class(el)[2], "] ", el$message)}
            )),
            collapse = ";\n"
          ),
          mist_statistics = ifelse(
            test = any(
              unique(sapply(
                X = has_issues, 
                FUN = function(el) {class(el)[2]}
              )) != "warning") ,
            yes = list(NA),
            no = list(mist_print(out_mist)$statistic) ## if only warning keep statistic
          ),
          mist_estimate = ifelse(
            test = any(
              unique(sapply(
                X = has_issues, 
                FUN = function(el) {class(el)[2]}
              )) != "warning") ,
            yes = list(NA),
            no = list(mist_print(out_mist)$estimate) ## if only warning keep estimates
          )
        ) 
      }
      stat_MiSTrare <- tidyr::unnest(data = out, cols = c(mist_statistics, mist_estimate)) 

    } else {
      ## analysis quanti
      
      ## Mist ## 
      has_issues <- try(tools::assertCondition(out_mist <- mist(
        y = y, 
        X = X, 
        G = G, 
        Z = Z, 
        model = "continuous"
      )), silent = TRUE)
      has_issues <- has_issues[
        sapply(has_issues, function(el) {length(intersect(class(el), c("warning", "error")))!=0})
      ]
      if (length(has_issues)==0) {
        out <- tibble::tibble(
          trait = Params["y_name"], 
          covariates = paste(covar_names, collapse = ";"),
          sample_size = nrow(Phenotypes),
          error = "None",
          mist_statistics = list(mist_print(out_mist)$statistic),
          mist_estimate = list(mist_print(out_mist)$estimate)
        ) 
      } else {
        out <- tibble::tibble(
          trait = Params["y_name"], 
          covariates = paste(covar_names, collapse = ";"),
          sample_size = nrow(Phenotypes),
          error = paste(
            unique(sapply(
              X = has_issues, 
              FUN = function(el) {paste0("[", class(el)[2], "] ", el$message)}
            )),
            collapse = ";\n"
          ),
          mist_statistics = ifelse(
            test = any(
              unique(sapply(
                X = has_issues, 
                FUN = function(el) {class(el)[2]}
              )) != "warning") ,
            yes = list(NA),
            no = list(mist_print(out_mist)$statistic)
          ), 
          mist_estimate = ifelse(
            test = any(
              unique(sapply(
                X = has_issues, 
                FUN = function(el) {class(el)[2]}
              )) != "warning") ,
            yes = list(NA),
            no =  list(mist_print(out_mist)$estimate)
          )
        )
      }
      stat_MiSTrare <- tidyr::unnest(data = out, cols = c(mist_statistics, mist_estimate)) 
      
      ## Add in clinicalData, mean Trait ~ rare variant count --
      
      tmp <- cbind.data.frame(Y = Phenotypes[, Params["y_name"]], G_rare)
      tmp$nb_mut = apply(X = tmp[,-1], MARGIN = 1, FUN = sum)
      
      meantrait_countmut <- tmp %>% 
        group_by(nb_mut) %>% 
        summarise(
          n = n(), 
          min_trait = min(Y),
          q25_trait = quantile(Y, probs = 0.25),
          median_trait = median(Y), 
          q75_trait = quantile(Y, probs = 0.75), 
          max_trait = max(Y),
          mean_trait = mean(Y), 
          sd_trait = sd(Y)
        )
      
      bx_plot <- ggplot(data = tmp, mapping = aes(x = as.factor(nb_mut), y = Y)) + 
        ggbeeswarm::geom_quasirandom(aes(color = as.factor(nb_mut)), width = 0.25, shape = 21) +
        geom_violin(fill = "transparent") +
        geom_boxplot(outlier.shape = NA, fill = "transparent", width = 0.25) + 
        theme(legend.position = "none", axis.ticks = element_blank()) + 
        scale_x_discrete(
          breaks = meantrait_countmut$nb_mut,
          labels = paste0(meantrait_countmut$nb_mut, "\n(N=", meantrait_countmut$n, ")")
        ) + 
        labs(
          title = paste0("Boxplot of ", Params["y_name"], " per number of carried risk alleles"),
          subtitle = paste0("MiST p.value.overall = ", round(x = stat_MiSTrare$p.value.overall, digits = 4)),
          x = paste0("Number of carried risk alleles\namong the cluster of ", ncol(G_rare)," rare variants"), 
          y = Params["y_name"], 
          caption = "N, the number of individuals carrying the allele(s)."
        ) 
      
      clinicalData <- append(x = clinicalData, values = list(bx_plot, meantrait_countmut))
    }
  
    stat_MiSTrare$nb_rare_var <- ncol(G_rare)  
  } else {
    message(paste("[MONO-GENE]", "After all QCs, not enough rare variants."))
    ## After all QC, not enough rare variants... Rare analyses impossible
    if (Params["binary"]) {
      stat_MiSTrare <- data.frame(matrix(rep(NA, 16), byrow = TRUE, nrow = 1))[-1, ]
      colnames(stat_MiSTrare) <- c(
        "trait", "covariates", "sample_size", "error", 
        "S.pi", "p.value.S.pi", "S.tau", "p.value.S.tau", "p.value.overall", 
        "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "OR", "nb_rare_var"
      )
    } else {
      stat_MiSTrare <- data.frame(matrix(rep(NA, 15), byrow = TRUE, nrow = 1))[-1, ]
      colnames(stat_MiSTrare) <- c(
        "trait", "covariates", "sample_size", "error",
        "S.pi", "p.value.S.pi", "S.tau", "p.value.S.tau", "p.value.overall", 
        "Pi_hat", "CI_2.5", "CI_97.5", "SE", "P_val", "nb_rare_var"
      )
    }
  }
  ## Output
  # stat_MiSTrare
  
  #### Effectif mutation  --------------------------------------------------------------------
  if (nrow(Phenotypes)>0) {
    if (Params["binary"]) {
      Effectif_mutation <- Effectif_mutation_YCC(
        Phenotypes = Phenotypes,
        G_mat = all_genotype_data,
        y_name = Params["y_name"]
      ) %>%
        select(-starts_with("genotype_Missing"))
      ## subset column missing_genotype because there is no samples with missing values in the final analyses
    } else {
      Effectif_mutation <- Effectif_mutation_Yquanti(
        G_mat = all_genotype_data[, -grep("IID", names(all_genotype_data)), drop = FALSE]
      )
    }
  } else {
    if (Params["binary"]) {
      Effectif_mutation <- structure(
        list(
          Variant = character(0), genotype_0_CTRL = character(0),
          genotype_1_CTRL = character(0), genotype_2_CTRL = character(0),
          genotype_0_CASE = character(0), genotype_1_CASE = character(0),
          genotype_2_CASE = character(0)
        ),
        row.names = integer(0), class = c(
          "tbl_df",
          "tbl", "data.frame"
        )
      )
    } else {
       Effectif_mutation <- structure(
        list(
          Variant = character(0), genotype_0 = character(0),
          genotype_1 = character(0), genotype_2 = character(0)
        ),
        row.names = integer(0), class = c(
          "tbl_df",
          "tbl", "data.frame"
        )
      )
    }
  }
  
  #### Output reshape ------------------------------------------------------------------------
  
  ## reshape QCobj
  if (QCdone) {
    xxqc <- f(Params) %>%
      as_tibble() %>%
      select(-GroupeGene)
    xxqc$QC1 <- list(f4(QC[[1]]))
    xxqc$QC2 <- list(f4(QC[[2]]))
    xxqc$QC3 <- list(f4(QC[[3]]))
    xxqc$QC4 <- list(f4(QC[[4]]))
    xxqc$QC5 <- list(f4(QC[[5]]))
    
  } else {
    xxqc <- f(Params) %>%
      as_tibble() %>%
      select(-GroupeGene)
  }
  
  if (nrow(Detail_Geno) == 0) {
    Annotation_qced <- left_join(
      x = Detail_Geno,
      y = all_quality_data %>%
        select(-Type_Variant) %>%
        slice(0),
      by = c("chr_pos")
    ) %>% 
      left_join(
        x = .,
        y = Effectif_mutation,
        by = "Variant"
      ) %>% 
      mutate(Samples_with_mutation = NA)
  } else {
    Annotation_qced <- left_join(
      x = Detail_Geno,
      y = all_quality_data,
      by = c("chr_pos", "Type_Variant")
    ) %>%
      left_join(
        x = .,
        y = Effectif_mutation,
        by = "Variant"
      )
  }
  
  res_rare <- stat_MiSTrare %>% 
    select(trait, covariates, sample_size, everything()) %>% 
    as.data.frame()
  
  ## Data analysed at the end
  if (nrow(Annotation_qced) > 0) {
    Annotation_qced <- Annotation_qced %>%
      mutate(
        Samples_with_mutation = map_chr(
          .x = Variant,
          .f = function(x) {
            all_genotype_data[, c(x, "IID"), drop = FALSE] %>%
              filter(get(x) != 0) %>%
              `[[`("IID") %>%
              paste(., collapse = ";")
          }
        )
      ) %>%
      select(Variant, Variant_Old, everything())
    ## For frequent variant Samples_with_mutation is too long for excel (limit of 32,767 characters)
    Annotation_qced[Annotation_qced$MAF > 0.05, "Samples_with_mutation"] <- "..."
    
    Annotation_qced <- Annotation_qced %>% 
      mutate(list_ind_mut = NULL, INFO = NULL)
  }
  
  #### SAVE IT ---------------------------------------------------------------------------------------
  
  if(!is.null(clinicalData) && length(clinicalData)>1){
    ## cas with quanti trait and ggplot to save
    ggsave(
      filename = paste0(
        OutDir, 
        Params["Gene_Symbol"], "/", 
        Params["y_name"], "_",  
        Params["analyse_LARGE_STRICT"], "_boxplot.png"
      ), 
      plot = clinicalData[2][[1]], 
      dpi = 300
    )
  
  }


  #### WRITE ALL RES to send by xlsx -----------------------------------------------------------------
  listf2 <- lapply(X = QC, FUN = f2)
  QC_table <- f3(listf2)
  
  if (nrow(Annotation_qced) > 0) {
    list_to_print <- list(
        "Annotation_qced" = Annotation_qced, 
        "QC_table" = QC_table, 
        "Bad_samples_tagged" = Bad_samples, 
        "clinicalData" = clinicalData[[1]]
      )
    if(!is.null(clinicalData) && length(clinicalData)>1){
      list_to_print$meantrait_countmut = clinicalData[3][[1]]
    }
    write_xlsx(
      x = list_to_print,
      path = paste0(
        OutDir, 
        Params["Gene_Symbol"], "/", 
        Params["y_name"], "_",  
        Params["analyse_LARGE_STRICT"], ".xlsx"
      ),
      col_names = TRUE
    )
    message(paste("[MONO-GENE]", "annotation_qced and QC_table ....  xlsx"))
  }
  
  ## to return in report ... 
  
  write_tsv(
    x = res_rare, 
    path = paste0(
      OutDir, 
      Params["Gene_Symbol"], "/", 
      Params["y_name"], "_",  
      Params["analyse_LARGE_STRICT"], "_resRare.tsv"
    ), 
    append = FALSE,
    col_names = TRUE
  )
  
  saveRDS(
    object = clinicalData, 
    file = paste0(
      OutDir, 
      Params["Gene_Symbol"], "/", 
      Params["y_name"], "_",  
      Params["analyse_LARGE_STRICT"], "_clinicalData.RDS"
    )
  )
  
  ## End of the analysis

} 

message(paste("[DATE-TIME]", Sys.time()))
