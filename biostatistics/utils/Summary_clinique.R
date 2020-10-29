Summary_clinique <- function(data, x, na_symbol = "NA") {
  continuous_trait <- function(x, digits = 3) {
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    if (is.na(x_mean)) {
      txt <- na_symbol
    } else {
      txt <- paste0(
        signif(x_mean, digits),
        " (", signif(x_sd, digits), ")"
      )
    }
    return(txt)
  }

  mean_sd_n <- function(x) {
    if (all(is.na(x))) {
      na_symbol ## "\\-"
    } else {
      paste(
        paste(
          format(mean(x, na.rm = TRUE), digits = 2, nsmall = 2, drop0trailing = FALSE),
          format(sd(x, na.rm = TRUE), digits = 2, nsmall = 2, drop0trailing = FALSE),
          sep = "±" # "\\$\\pm\\$"
        ),
        paste0(
          "(n=",
          sum(!is.na(x)),
          ")"
        )
      )
    }
  }
  ## drop0trailing = FALSE keep 0 in right

  out <- data %>%
    mutate(Status = get(x)) %>%
    filter(!is.na(Status)) %>%
    group_by(Status) %>%
    dplyr::summarise(
      N = format(dplyr::n(), big.mark = ",", scientific = FALSE),
      AGE = mean_sd_n(AGE),
      BMI = mean_sd_n(BMI),
      SEX = paste0(
        "M:", format(sum(SEX == 1), big.mark = ",", scientific = FALSE),
        " / F:", format(sum(SEX == 2), big.mark = ",", scientific = FALSE)
      ),
      FG = mean_sd_n(FG)
    ) %>%
    mutate(Status = c("0" = "Controls", "1" = "Cases")[as.character(Status)]) %>%
    as.data.frame() %>%
    column_to_rownames(var = "Status") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = x)

  if (ncol(out) == 1) {
    # out <- out %>% mutate(Data = na_symbol)
    out <- NULL
  }

  return(out)
}
