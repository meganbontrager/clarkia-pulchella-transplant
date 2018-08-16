# do edge populations have an advantage once climate is accounted for?

library(tidyverse)
library(cowplot)
library(glmmTMB)
library(ggeffects)


# load data ----------------------------------------------------

# load data and filter to within-population crosses
plants_wi =  read.csv("data/all_data.csv") %>% filter(type == "WI")

# create local/foreign column where focal pops in their home sites are local
plants_wi$local_foreign = ifelse((plants_wi$sirepop == "AQ" & plants_wi$site == "AQ")|(plants_wi$sirepop == "AD" & plants_wi$site == "AD"), "local", "foreign")

plants_wi_seeds = filter(plants_wi, !is.na(total_est_seeds))


# lifetime fitness -----------------------------------------------

seedsall.lf.mod = glmmTMB(total_est_seeds ~ local_foreign + abs_tave_diff_sep_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled
                          + (1|site/blk) + (1|sirepop), 
                          ziformula = ~ local_foreign + abs_tave_diff_sep_jul_scaled,
                          data = plants_wi_seeds, 
                          family = nbinom2)
summary(seedsall.lf.mod) # no significant difference between locals and foreigns


# other lifestages with size as covariate, truncated datasets -----------------------------------------------------

plants_wi$nov_size_scaled = as.vector(scale(plants_wi$nov_size))
plants_wi$mar_size_scaled = as.vector(scale(plants_wi$mar_size))

germ.lf.mod = glmmTMB(nov_germ ~ local_foreign + abs_tave_diff_sep_nov_scaled
                      + (1|site/blk) + (1|sirepop), 
                      family = binomial, 
                      data = plants_wi)
summary(germ.lf.mod)

novsize.lf.mod = glmmTMB(nov_size ~ local_foreign + abs_tave_diff_sep_nov_scaled
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian, 
                         data = filter(plants_wi, !is.na(nov_size_scaled)))
summary(novsize.lf.mod)

marsurv.lf.mod = glmmTMB(mar_surv ~ local_foreign + nov_size_scaled + abs_tave_diff_dec_mar_scaled
                         + (1|sirepop) + (1|site), 
                         family = binomial, 
                         data = filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.lf.mod)

marsize.lf.mod = glmmTMB(mar_size ~ local_foreign + nov_size_scaled + abs_tave_diff_dec_mar_scaled
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian,
                         data = filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.lf.mod)

fruitcount.lf.mod = glmmTMB(fruit_count ~ local_foreign + mar_size_scaled + abs_tave_diff_apr_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled
                            + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                            family = nbinom2, 
                            ziformula = ~ 1,
                            data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.lf.mod)

seeds.lf.mod = glmmTMB(total_est_seeds ~ local_foreign + mar_size_scaled + abs_tave_diff_apr_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled
                       + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                       family = nbinom2, 
                       ziformula = ~ 1,
                       data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.lf.mod)




