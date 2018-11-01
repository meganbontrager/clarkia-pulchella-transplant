# Does gene flow help or hurt edge populations?
# referenced in Figure 4, Table S4

# libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(GGally)
library(glmmTMB) # version 0.2.2.0
library(ggeffects) # version 0.3.0



# load data ---------------------------------------------------------------

# load full data
plants = read.csv("data/Bontrager_transplant_data.csv") 

# generate local/foreign column
plants$local_foreign = ifelse((plants$sirepop == "AQ" & plants$site == "AQ")|(plants$sirepop == "AD" & plants$site == "AD"), "local", "foreign")

# filter to plants with reproduction estimates
plants_seeds = plants %>% filter(!is.na(total_est_seeds))



# overall fitness ---------------------------------------------------------

seedsall.gf.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + type + abs_ppt_mm_diff_midparent_apr_jul_scaled
                          + (1|site/blk) + (1|sirepop), 
                          ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + type,
                          data = plants_seeds, 
                          family = nbinom2)
summary(seedsall.gf.mod) 
nrow(plants_seeds)
# 27312


seedsall.gf.noppt.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + type
                                + (1|site/blk) + (1|sirepop), 
                                ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + type,
                                data = plants_seeds, 
                                family = nbinom2)
summary(seedsall.gf.noppt.mod) 



seedsall.gf.notype.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + 
                                 abs_ppt_mm_diff_midparent_apr_jul_scaled
                                 + (1|site/blk) + (1|sirepop), 
                                 ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
                                 data = plants_seeds, 
                                 family = nbinom2)
summary(seedsall.gf.notype.mod) 



# other lifestages  -----------------------------------------------------

# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

# generate scaled size covariate
plants$nov_size_scaled = as.vector(scale(plants$nov_size))
plants$mar_size_scaled = as.vector(scale(plants$mar_size))

germ.gf.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled + type
                      + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = binomial, 
                      data = filter(plants, !is.na(nov_germ)))
summary(germ.gf.mod)
nrow(filter(plants, !is.na(nov_germ)))
# 32361

novsize.gf.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled + type
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = gaussian,
                      data = filter(plants, nov_germ == 1, !is.na(nov_size)))
summary(novsize.gf.mod)
nrow(filter(plants, nov_germ == 1, !is.na(nov_size)))
# 8605

marsurv.gf.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + type + nov_size_scaled
                         + (1|site) + (1|sirepop), 
                       family = binomial, 
                       data = filter(plants, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.gf.mod)
nrow(filter(plants, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
# 8531

marsize.gf.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + type + nov_size_scaled
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = gaussian,
                      data = filter(plants, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.gf.mod)
nrow(filter(plants, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
# 5249

fruitcount.gf.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + type + mar_size_scaled
                            + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.gf.mod)
nrow(filter(plants, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
# 6231

seeds.gf.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + type + mar_size_scaled
                       + (1|site/blk) + (1|sirepop), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.gf.mod)
nrow(filter(plants, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
# 6231


# build a results table -----------------------------------------

mod_list = ls(pattern = "*.mod")

for (i in 1:length(mod_list)){
  mod = get(as.character(mod_list[i]))
  mod_sum = summary(get(as.character(mod_list[i])))
  c = rownames_to_column(data.frame(t(mod_sum$coefficients$cond)))
  c1 = cbind(model = paste(mod_list[i]), part = "cond", c) 
  if (is.null(mod_sum$coefficients$zi)) {
    z = NA
    z1 = NA
    results = c1
  } else {
    z = rownames_to_column(data.frame(t(mod_sum$coefficients$zi)))
    z1 = cbind(model = paste(mod_list[i]), part = "zi", z) 
    results = bind_rows(c1, z1)
  }
  if (i == 1) results_new = results else results_new = bind_rows(results_old, results)
  results_old = results_new
}

results_gf = results_new %>% 
  unite(model, model, part, remove = TRUE, sep = ".") %>% 
  gather(var, value, c(3:11)) %>% 
  mutate(var = ifelse(var %in% c("abs_tave_diff_midparent_apr_jul_scaled", "abs_tave_diff_midparent_dec_mar_scaled", 
                                 "abs_tave_diff_midparent_sep_jul_scaled", "abs_tave_diff_midparent_sep_nov_scaled", 
                                 "abs_tave_diff_midparent_apr_may_scaled"), "temp", var)) %>%
  mutate(var = ifelse(var %in% c("abs_ppt_mm_diff_midparent_apr_jul_scaled", "abs_ppt_mm_diff_midparent_apr_may_scaled"), "ppt", var)) %>% 
  filter(!is.na(value)) %>%
  unite(variable, var, rowname, sep = ":") %>% 
  spread(variable, value)

results_round = rapply(object = results_gf, f = round, classes = "numeric", how = "replace", digits = 3) 

# write.csv(results_round, "results/gf_table.csv", row.names = FALSE)



