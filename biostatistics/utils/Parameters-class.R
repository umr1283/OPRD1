######################################################################
####################### Class OPRD1_Parameters #######################
############################## Creation ##############################
######################################################################


### Class definition ###
setClass(
  Class = "OPRD1_Parameters", 
  representation = representation(
    FolderOut = "character", 
    GroupeGene = "character", 
    Gene_Symbol = "character", 
    ENST = "character", 
    analyse_LARGE_STRICT = "character", 
    threshold = "numeric", 
    ExludeSynonymous = "logical", 
    cluster = "character", 
    y_name = "character", 
    binary = "logical", 
    SNP_INDEL = "character"
  ), 
  prototype = prototype(
    FolderOut = character(), 
    GroupeGene = character(), 
    Gene_Symbol = character(), 
    ENST = character(), 
    analyse_LARGE_STRICT = character(), 
    threshold = numeric(), 
    ExludeSynonymous = logical(), 
    cluster = character(), 
    y_name = character(), 
    binary = logical(), 
    SNP_INDEL = character()
  )# , 
  # validity = function (object) {
    # cat("**** validity OPRD1_Parameters <empty> ****\n")
    # return(TRUE)
  # }
)


### Constructor ###
setGeneric(name = "new.OPRD1_Parameters", def = function (FolderOut, GroupeGene, Gene_Symbol, ENST, analyse_LARGE_STRICT, threshold, ExludeSynonymous, cluster, y_name, binary, SNP_INDEL) {standardGeneric("new.OPRD1_Parameters")})
setMethod(f = "new.OPRD1_Parameters", signature = c("missing", "missing", "missing", "missing", "missing", "missing", "missing", "missing", "missing", "missing", "missing"), definition = function (FolderOut, GroupeGene, Gene_Symbol, ENST, analyse_LARGE_STRICT, threshold, ExludeSynonymous, cluster, y_name, binary, SNP_INDEL) {new("OPRD1_Parameters")})
setMethod(f = "new.OPRD1_Parameters", signature = c("ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY"), definition = function (FolderOut, GroupeGene, Gene_Symbol, ENST, analyse_LARGE_STRICT, threshold, ExludeSynonymous, cluster, y_name, binary, SNP_INDEL) {
  if (missing(FolderOut)) {FolderOut <- character()} else {}
  if (missing(GroupeGene)) {GroupeGene <- character()} else {}
  if (missing(Gene_Symbol)) {Gene_Symbol <- character()} else {}
  if (missing(ENST)) {ENST <- character()} else {}
  if (missing(analyse_LARGE_STRICT)) {analyse_LARGE_STRICT <- character()} else {}
  if (missing(threshold)) {threshold <- numeric()} else {}
  if (missing(ExludeSynonymous)) {ExludeSynonymous <- logical()} else {}
  if (missing(cluster)) {cluster <- character()} else {}
  if (missing(y_name)) {y_name <- character()} else {}
  if (missing(binary)) {binary <- logical()} else {}
  if (missing(SNP_INDEL)) {SNP_INDEL <- character()} else {}
  return(new("OPRD1_Parameters", FolderOut = FolderOut, GroupeGene = GroupeGene, Gene_Symbol = Gene_Symbol, ENST = ENST, analyse_LARGE_STRICT = analyse_LARGE_STRICT, threshold = threshold, ExludeSynonymous = ExludeSynonymous, cluster = cluster, y_name = y_name, binary = binary, SNP_INDEL = SNP_INDEL))
})


### Is ###
setGeneric(name = "is.OPRD1_Parameters", def = function (object) {standardGeneric("is.OPRD1_Parameters")})
setMethod(f = "is.OPRD1_Parameters", signature = "ANY", definition = function (object) {
  if (length(object)>1) {
    return(sapply(object, is.OPRD1_Parameters))
  } else {
    if (class(object) == "OPRD1_Parameters") {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
})


### Show ###
setMethod(f = "show", signature = "OPRD1_Parameters", definition = function (object){
  showSlot <- function (slot) {
    sNames <- gsub("^[^@]*@(.*)", "\\1", slot)
    eSlot <- eval(parse(text = slot))
    tmp <- switch(EXPR = class(eSlot),
      "matrix" = {
        cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
        if (all(dim(eSlot)==0)) {
          cat("NA")
        } else {
          cat("\n")
          nrowShow <- seq(min(5, nrow(eSlot)))
          ncolShow <- seq(min(5, ncol(eSlot)))
          shortObject <- eSlot[nrowShow, ncolShow]
          if (is.null(rownames(shortObject))) {
            rownames(shortObject) <- seq(nrow(shortObject))
          } else {}
          if (is.null(colnames(shortObject))) {
            colnames(shortObject) <- seq(ncol(shortObject))
          } else {}
          resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
          if (nrow(shortObject)!=nrow(eSlot)) {
            resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function (iCol) {paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")})))
          } else {}
          if (ncol(shortObject)!=ncol(eSlot)) {
            resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat)-1)))
          } else {}
          cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
        }
        cat("\n")
      },
      "data.frame" = {
        cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
        if (all(dim(eSlot)==0)) {
          cat(" NA")
        } else {
          cat("\n")
          nrowShow <- seq(min(5, nrow(eSlot)))
          ncolShow <- seq(min(5, ncol(eSlot)))
          shortObject <- eSlot[nrowShow, ncolShow]
          if (is.null(rownames(shortObject))) {
            rownames(shortObject) <- seq(nrow(shortObject))
          } else {}
          if (is.null(colnames(shortObject))) {
            colnames(shortObject) <- seq(ncol(shortObject))
          } else {}
          resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
          if (nrow(shortObject)!=nrow(eSlot)) {
            resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function (iCol) {paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")})))
          } else {}
          if (ncol(shortObject)!=ncol(eSlot)) {
            resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat)-1)))
          } else {}
          cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
        }
        cat("\n")
      },
      "numeric" = {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot)>1) {
            cat(paste0("[", length(eSlot), "] ", paste0(format(head(eSlot), digits = 4), collapse = " ")))
          } else {
            cat(format(eSlot, digits = 4))
          }
        }
        cat("\n")
      },
      "character" = {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot)>1) {
            cat("[", length(eSlot), "] \"", paste0(head(eSlot), collapse = "\" \""), "\"", sep = "")
          } else {
            cat(paste0("\"", eSlot, "\""))
          }
        }
        cat("\n")
      },
      {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot)>1) {
            cat(paste0("[", length(eSlot), "] ", paste0(head(eSlot), collapse = " ")))
          } else {
            cat(eSlot)
          }
        }
        cat("\n")
      }
    )
    return(invisible())
  }
  showObject <- function (object) {
    cat("  ~~~ Class:", class(object), "~~~\n")
    sNames <- paste0("object@", slotNames(object))
    trash <- sapply(sNames, showSlot)
    return(invisible())
  }
  showObject(object)
  return(invisible(object))
})


### Getteur ###
setMethod(f = "[", signature = "OPRD1_Parameters", definition = function (x, i, j, drop){
  switch(EXPR = i, 
    "FolderOut" = {
      if (missing(j)) {
        return(x@FolderOut)
      } else {
        if (j>length(x@FolderOut)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@FolderOut[j])
        }
      }
    }, 
    "GroupeGene" = {
      if (missing(j)) {
        return(x@GroupeGene)
      } else {
        if (j>length(x@GroupeGene)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@GroupeGene[j])
        }
      }
    }, 
    "Gene_Symbol" = {
      if (missing(j)) {
        return(x@Gene_Symbol)
      } else {
        if (j>length(x@Gene_Symbol)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@Gene_Symbol[j])
        }
      }
    }, 
    "ENST" = {
      if (missing(j)) {
        return(x@ENST)
      } else {
        if (j>length(x@ENST)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@ENST[j])
        }
      }
    }, 
    "analyse_LARGE_STRICT" = {
      if (missing(j)) {
        return(x@analyse_LARGE_STRICT)
      } else {
        if (j>length(x@analyse_LARGE_STRICT)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@analyse_LARGE_STRICT[j])
        }
      }
    }, 
    "threshold" = {
      if (missing(j)) {
        return(x@threshold)
      } else {
        if (j>length(x@threshold)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@threshold[j])
        }
      }
    }, 
    "ExludeSynonymous" = {
      if (missing(j)) {
        return(x@ExludeSynonymous)
      } else {
        if (j>length(x@ExludeSynonymous)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@ExludeSynonymous[j])
        }
      }
    }, 
    "cluster" = {
      if (missing(j)) {
        return(x@cluster)
      } else {
        if (j>length(x@cluster)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@cluster[j])
        }
      }
    }, 
    "y_name" = {
      if (missing(j)) {
        return(x@y_name)
      } else {
        if (j>length(x@y_name)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@y_name[j])
        }
      }
    }, 
    "binary" = {
      if (missing(j)) {
        return(x@binary)
      } else {
        if (j>length(x@binary)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@binary[j])
        }
      }
    }, 
    "SNP_INDEL" = {
      if (missing(j)) {
        return(x@SNP_INDEL)
      } else {
        if (j>length(x@SNP_INDEL)) {
          stop("[OPRD1_Parameters:get] indice out of limits")
        } else {
          return(x@SNP_INDEL[j])
        }
      }
    }, 
    stop("[OPRD1_Parameters:get] ", i, " is not a \"OPRD1_Parameters\" slot")
  )
})


### Setteur ###
setMethod(f = "[<-", signature = "OPRD1_Parameters", definition = function (x, i, j, value){
    switch(EXPR = i, 
    "FolderOut" = {
      if (missing(j)) {
        x@FolderOut <- value
      } else {
        if (j>length(x@FolderOut)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@FolderOut[j] <- value
        }
      }
    }, 
    "GroupeGene" = {
      if (missing(j)) {
        x@GroupeGene <- value
      } else {
        if (j>length(x@GroupeGene)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@GroupeGene[j] <- value
        }
      }
    }, 
    "Gene_Symbol" = {
      if (missing(j)) {
        x@Gene_Symbol <- value
      } else {
        if (j>length(x@Gene_Symbol)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@Gene_Symbol[j] <- value
        }
      }
    }, 
    "ENST" = {
      if (missing(j)) {
        x@ENST <- value
      } else {
        if (j>length(x@ENST)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@ENST[j] <- value
        }
      }
    }, 
    "analyse_LARGE_STRICT" = {
      if (missing(j)) {
        x@analyse_LARGE_STRICT <- value
      } else {
        if (j>length(x@analyse_LARGE_STRICT)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@analyse_LARGE_STRICT[j] <- value
        }
      }
    }, 
    "threshold" = {
      if (missing(j)) {
        x@threshold <- value
      } else {
        if (j>length(x@threshold)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@threshold[j] <- value
        }
      }
    }, 
    "ExludeSynonymous" = {
      if (missing(j)) {
        x@ExludeSynonymous <- value
      } else {
        if (j>length(x@ExludeSynonymous)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@ExludeSynonymous[j] <- value
        }
      }
    }, 
    "cluster" = {
      if (missing(j)) {
        x@cluster <- value
      } else {
        if (j>length(x@cluster)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@cluster[j] <- value
        }
      }
    }, 
    "y_name" = {
      if (missing(j)) {
        x@y_name <- value
      } else {
        if (j>length(x@y_name)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@y_name[j] <- value
        }
      }
    }, 
    "binary" = {
      if (missing(j)) {
        x@binary <- value
      } else {
        if (j>length(x@binary)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@binary[j] <- value
        }
      }
    }, 
    "SNP_INDEL" = {
      if (missing(j)) {
        x@SNP_INDEL <- value
      } else {
        if (j>length(x@SNP_INDEL)) {
          stop("[OPRD1_Parameters:set] indice out of limits")
        } else {
          x@SNP_INDEL[j] <- value
        }
      }
    }, 
    stop("[OPRD1_Parameters:set] ", i, " is not a \"OPRD1_Parameters\" slot")
  )
  validObject(x)
  return(invisible(x))
})


### Summary ###
setMethod(f = "summary", signature = "OPRD1_Parameters", definition = function (object){
  if (missing(object)){
    stop("[OPRD1_Parameters:summary] \"object\" is missing", call. = FALSE)
    return(invisible())
  } else {}
  warning("[OPRD1_Parameters:summary] No summary method defined for \"OPRD1_Parameters\" object!", call. = FALSE)
  return(invisible(object))
})
