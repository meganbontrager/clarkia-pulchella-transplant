# categorical comparison of local vs. foreign warmer vs. foreign cooler
# referenced in supplemental analyses 1, figure S5


# libraries ---------------------------------------------------------------

library(tidyverse)
library(glmmTMB) # version 0.2.2.0


# data --------------------------------------------------------------------

plants = read.csv("data/Bontrager_transplant_data.csv") %>% filter(type == "WI")

# generate local/foreign warmer/foreign cooler column
plants$local_foreign_warmer = as.factor(ifelse((plants$sirepop == "AQ" & plants$site == "AQ")|(plants$sirepop == "AD" & plants$site == "AD"), "local", 
                                               ifelse((plants$site == "AQ" & plants$tave_5180_sep_jul > 5.13)|(plants$site == "AD" & plants$tave_5180_sep_jul > 4.50), "foreign_warmer", "foreign_cooler")))

plants$local_foreign_warmer = factor(plants$local_foreign_warmer, levels = c("local", "foreign_warmer", "foreign_cooler"))

table(plants$local_foreign_warmer)

plants_seeds = filter(plants, !is.na(total_est_seeds))

table(plants_seeds$local_foreign_warmer)

seedsall.pops.mod = glmmTMB(total_est_seeds ~ local_foreign_warmer
                            + (1|site/blk) + (1|sirepop),
                            ziformula = ~ local_foreign_warmer,
                            data = plants_seeds, 
                            family = nbinom2)
summary(seedsall.pops.mod)
