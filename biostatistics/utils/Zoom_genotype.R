Zoom_Genotype <- function(MatriceG) {
  if (ncol(MatriceG) > 0) {
    if (class(MatriceG) != "data.frame") {
      MatriceG <- as.data.frame(MatriceG)
    }

    computeMAF <- function(x, coding = c(-1, 0, 1, 2)) {
      if (all(x == -1)) {
        maf <- 0
      } else {
        dta.eff <- table(factor(x, levels = coding))[c("0", "1", "2")]
        maf <- sum(dta.eff * c(0, 1, 2)) / sum(dta.eff * 2)
      }
      return(maf)
    }
    
    
    nb_cores <- min(20, ncol(MatriceG))
    # set_ncores(nb_cores) ## message
    Mat_eff <- mclapply(
      MatriceG,
      mc.cores = nb_cores,
      mc.preschedule = FALSE,
      function(icol) {
        if (computeMAF(icol) > 0.5) {
          icol <- abs(icol - 2)
          icol[icol == 3] <- -1
          dta.properref <- FALSE
          # Reverse REF and ALT to have Ref as Maj_allele and Alt as Minor_Allele
        } else {
          dta.properref <- TRUE
        }

        dta.eff <- table(factor(icol, levels = c(-1, 0, 1, 2)))
        names(dta.eff) <- paste("G", names(dta.eff), sep = "_")
        names(dta.eff) <- gsub("G_-1", "Missing_genotype", names(dta.eff))
        Mat_eff <- c(dta.eff, MAF = computeMAF(icol), ProperAlleleREF = dta.properref)
        return(Mat_eff)
      }
    )
    # close_ncores(nb_cores) ## cores free

    Mat_eff <- as.data.frame(do.call("rbind", Mat_eff))
    Mat_eff[, "Variant"] <- rownames(Mat_eff)
    Mat_eff <- Mat_eff[, c(ncol(Mat_eff), seq_len(ncol(Mat_eff) - 1))] ## reordering column
    rownames(Mat_eff) <- NULL
    Mat_eff$VarClean <- ifelse(Mat_eff$ProperAlleleREF == 1, Mat_eff$Variant, gsub("([0-9]*_[0-9]*_)(.*)_(.*)", "\\1\\3_\\2", Mat_eff$Variant))
    ## check Mat_eff[Mat_eff$ProperAlleleREF==0,]

    Mat_eff$chr_pos <- gsub("([0-9]*_[0-9]*)_(.*)_(.*)", "\\1", Mat_eff$Variant)
  } else {
    Mat_eff <- data.frame(matrix(nrow = 0, ncol = 9))
    names(Mat_eff) <- c("Variant", "Missing_genotype", "G_0", "G_1", "G_2", "MAF", "ProperAlleleREF", "VarClean", "chr_pos")
  }
  return(Mat_eff)
}
