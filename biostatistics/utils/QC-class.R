######################################################################
########################### Class OPRD1_QC ###########################
############################## Creation ##############################
######################################################################


### Class definition ###
setClass(
  Class = "OPRD1_QC", 
  representation = representation(
    Step = "numeric", 
    Why = "character", 
    N_origin = "numeric", 
    Excluded = "character", 
    N_excluded = "numeric", 
    Type = "character", 
    Comment = "character"
  ), 
  prototype = prototype(
    Step = numeric(), 
    Why = character(), 
    N_origin = numeric(), 
    Excluded = character(), 
    N_excluded = numeric(), 
    Type = character(), 
    Comment = character()
  )# , 
  # validity = function (object) {
    # cat("**** validity OPRD1_QC <empty> ****\n")
    # return(TRUE)
  # }
)


### Constructor ###
setGeneric(name = "new.OPRD1_QC", def = function (Step, Why, N_origin, Excluded, N_excluded, Type, Comment) {standardGeneric("new.OPRD1_QC")})
setMethod(f = "new.OPRD1_QC", signature = c("missing", "missing", "missing", "missing", "missing", "missing", "missing"), definition = function (Step, Why, N_origin, Excluded, N_excluded, Type, Comment) {new("OPRD1_QC")})
setMethod(f = "new.OPRD1_QC", signature = c("ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY"), definition = function (Step, Why, N_origin, Excluded, N_excluded, Type, Comment) {
  if (missing(Step)) {Step <- numeric()} else {}
  if (missing(Why)) {Why <- character()} else {}
  if (missing(N_origin)) {N_origin <- numeric()} else {}
  if (missing(Excluded)) {Excluded <- character()} else {}
  if (missing(N_excluded)) {N_excluded <- numeric()} else {}
  if (missing(Type)) {Type <- character()} else {}
  if (missing(Comment)) {Comment <- character()} else {}
  return(new("OPRD1_QC", Step = Step, Why = Why, N_origin = N_origin, Excluded = Excluded, N_excluded = N_excluded, Type = Type, Comment = Comment))
})


### Is ###
setGeneric(name = "is.OPRD1_QC", def = function (object) {standardGeneric("is.OPRD1_QC")})
setMethod(f = "is.OPRD1_QC", signature = "ANY", definition = function (object) {
  if (length(object)>1) {
    return(sapply(object, is.OPRD1_QC))
  } else {
    if (class(object) == "OPRD1_QC") {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
})


### Show ###
setMethod(f = "show", signature = "OPRD1_QC", definition = function (object){
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
setMethod(f = "[", signature = "OPRD1_QC", definition = function (x, i, j, drop){
  switch(EXPR = i, 
    "Step" = {
      if (missing(j)) {
        return(x@Step)
      } else {
        if (j>length(x@Step)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@Step[j])
        }
      }
    }, 
    "Why" = {
      if (missing(j)) {
        return(x@Why)
      } else {
        if (j>length(x@Why)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@Why[j])
        }
      }
    }, 
    "N_origin" = {
      if (missing(j)) {
        return(x@N_origin)
      } else {
        if (j>length(x@N_origin)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@N_origin[j])
        }
      }
    }, 
    "Excluded" = {
      if (missing(j)) {
        return(x@Excluded)
      } else {
        if (j>length(x@Excluded)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@Excluded[j])
        }
      }
    }, 
    "N_excluded" = {
      if (missing(j)) {
        return(x@N_excluded)
      } else {
        if (j>length(x@N_excluded)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@N_excluded[j])
        }
      }
    }, 
    "Type" = {
      if (missing(j)) {
        return(x@Type)
      } else {
        if (j>length(x@Type)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@Type[j])
        }
      }
    }, 
    "Comment" = {
      if (missing(j)) {
        return(x@Comment)
      } else {
        if (j>length(x@Comment)) {
          stop("[OPRD1_QC:get] indice out of limits")
        } else {
          return(x@Comment[j])
        }
      }
    }, 
    stop("[OPRD1_QC:get] ", i, " is not a \"OPRD1_QC\" slot")
  )
})


### Setteur ###
setMethod(f = "[<-", signature = "OPRD1_QC", definition = function (x, i, j, value){
    switch(EXPR = i, 
    "Step" = {
      if (missing(j)) {
        x@Step <- value
      } else {
        if (j>length(x@Step)) {
          stop("[OPRD1_QC:set] indice out of limits")
        } else {
          x@Step[j] <- value
        }
      }
    }, 
    "Why" = {
      if (missing(j)) {
        x@Why <- value
      } else {
        if (j>length(x@Why)) {
          stop("[OPRD1_QC:set] indice out of limits")
        } else {
          x@Why[j] <- value
        }
      }
    }, 
    "N_origin" = {
      if (missing(j)) {
        x@N_origin <- value
      } else {
        if (j>length(x@N_origin)) {
          stop("[OPRD1_QC:set] indice out of limits")
        } else {
          x@N_origin[j] <- value
        }
      }
    }, 
    "Excluded" = {
      if (missing(j)) {
        x@Excluded <- value
      } else {
        if (j>length(x@Excluded)) {
          stop("[OPRD1_QC:set] indice out of limits")
        } else {
          x@Excluded[j] <- value
        }
      }
    }, 
    "N_excluded" = {
      if (missing(j)) {
        x@N_excluded <- value
      } else {
        if (j>length(x@N_excluded)) {
          stop("[OPRD1_QC:set] indice out of limits")
        } else {
          x@N_excluded[j] <- value
        }
      }
    }, 
    "Type" = {
      if (missing(j)) {
        x@Type <- value
      } else {
        if (j>length(x@Type)) {
          stop("[OPRD1_GENES_QC:set] indice out of limits")
        } else {
          x@Type[j] <- value
        }
      }
    }, 
    "Comment" = {
      if (missing(j)) {
        x@Comment <- value
      } else {
        if (j>length(x@Comment)) {
          stop("[OPRD1_GENES_QC:set] indice out of limits")
        } else {
          x@Comment[j] <- value
        }
      }
    }, 
    stop("[OPRD1_GENES_QC:set] ", i, " is not a \"OPRD1_GENES_QC\" slot")
  )
  validObject(x)
  return(invisible(x))
})


### Summary ###
setMethod(f = "summary", signature = "OPRD1_GENES_QC", definition = function (object){
  if (missing(object)){
    stop("[OPRD1_GENES_QC:summary] \"object\" is missing", call. = FALSE)
    return(invisible())
  } else {}
  warning("[RADiO_QC:summary] No summary method defined for \"RADiO_QC\" object!", call. = FALSE)
  return(invisible(object))
})
