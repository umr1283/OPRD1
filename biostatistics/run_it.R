#!/usr/bin/env Rscript

## ex : 
## nohup Rscript OPRD1_multiple_analyses.R GPCR GOF CC_OBESITE TRUE & 
## nohup Rscript OPRD1_multiple_analyses.R GPCR LOF CC_OBESITE TRUE & 
## nohup Rscript OPRD1_multiple_analyses.R GPCR ALL CC_OBESITE TRUE & 
##... 

options(stringsAsFactors = FALSE)
library(readxl)
library(tidyverse)

#### List it ####

Analyses_plan <- read_excel(
  path = "/disks/PROJECT/OPRD1/Docs/Analyses_plan_OPRD1.xlsx", 
  col_names = TRUE,
  na = "NA",
)

Add_LOF_modify <- Analyses_plan %>%
  ## add Version LOFgene+LOFfunc 
  filter(Version %in% "LOF") %>% 
  mutate(
    Version = "LOFgene+LOFfunc", 
    Model = map2_chr(.x = Model, .y = Version, .f = ~gsub(pattern = "LOF", replacement = .y, x = .x))
  ) 

Analyses_plan <- bind_rows(Analyses_plan, Add_LOF_modify) %>% 
  mutate(
    Binaire = as.logical(Binaire), 
    cmd = paste(
      "Rscript OPRD1_multiple_analyses.R GPCR ", Version, Y, Binaire, "& ")
  ) %>% 
  arrange(Y)

# writexl::write_xlsx(
#   x = Analyses_plan, 
#   path = "/disks/PROJECT/OPRD1/Docs/Analyses_plan_OPRD1.xlsx", 
#   col_names = TRUE
# )

# file.copy(
#   from = "/disks/PROJECT/OPRD1/Docs/Analyses_plan_OPRD1.xlsx", 
#   to = "/disks/PROJECT/OPRD1/Data/OPRD1_multiple_analyses/xlsx_details/Analyses_plan_OPRD1.xlsx", 
#   overwrite = TRUE
# )


## if want to run only the analyses "LOFgene+LOFfunc" ## example : 
# Analyses_plan <- Analyses_plan %>% filter(Version %in% "LOFgene+LOFfunc")
# Analyses_plan <- Analyses_plan %>% filter(!Version %in% "LOFgene+LOFfunc")


#### run it ####
walk(.x = Analyses_plan$cmd, .f = ~system(command = .x, intern = TRUE))

#### Make the archive ####
# "tar zcvf archive.tar.gz archive/"
