# Do the effects of gene flow depend upon the genetic differentiation between focal and donor populations?
# referenced in figure 5, table S5, table S6

library(tidyverse)
library(glmmTMB) # version 0.2.2.0


# load data ---------------------------------------------------------------

# load and filter to between-population crosses
plants_gf = read.csv("data/Bontrager_transplant_data.csv") %>% filter(type == "GF")
nrow(plants_gf)
# 16236

# filter to rows with seeds
plants_gf_seeds = plants_gf %>% filter(!is.na(total_est_seeds))
nrow(plants_gf_seeds)
# 13627



# correlation of predictors -----------------------------------------------

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_sep_jul_scaled", "fst_wc_scaled", "abs_ppt_mm_diff_midparent_apr_jul_scaled")])
ggpairs(preds_temp)
# ppt is correlated with fst at 0.635

preds_temp = unique(plants_gf[,c("abs_ppt_mm_diff_midparent_apr_jul_scaled", "fst_wc_scaled")])
ggpairs(preds_temp)

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_sep_nov_scaled", "fst_wc_scaled")])
ggpairs(preds_temp)

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_dec_mar_scaled", "fst_wc_scaled")])
ggpairs(preds_temp)

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_apr_jul_scaled", "fst_wc_scaled", "abs_ppt_mm_diff_midparent_apr_jul_scaled")])
ggpairs(preds_temp)


# lifetime fitness --------------------------------------------------------

# full model
seedsall.fst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + 
                          abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled
                                    + (1|site/blk) + (1|sirepop), 
                                    ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + fst_wc_scaled,
                                    data = plants_gf_seeds, 
                                    family = nbinom2)
summary(seedsall.fst.mod) 
# fst and temp significant in zero-inflation part

# univariate models
seedsall.fst.fstonly.mod = glmmTMB(total_est_seeds ~ fst_wc_scaled
                                   + (1|site/blk) + (1|sirepop), 
                                   ziformula =  ~ fst_wc_scaled ,
                                   data = plants_gf_seeds, 
                                   family = nbinom2)
summary(seedsall.fst.fstonly.mod) 
# fst significant in cond. and zi parts

seedsall.fst.temponly.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled
                                    + (1|site/blk) + (1|sirepop), 
                                    ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
                                    data = plants_gf_seeds, 
                                    family = nbinom2)
summary(seedsall.fst.temponly.mod) 
# temp significant in zi part

seedsall.fst.preciponly.mod = glmmTMB(total_est_seeds ~ abs_ppt_mm_diff_midparent_apr_jul_scaled
                                      + (1|site/blk) + (1|sirepop), 
                                      ziformula =  ~ 1,
                                      data = plants_gf_seeds, 
                                      family = nbinom2)
summary(seedsall.fst.preciponly.mod) 
# precip not significant 

# precip is not significant in univariate models (below) and is co-linear with fst (high fst == high difference in precip), so trying without it.
seedsall.fst.noprecip.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + fst_wc_scaled
+ (1|site/blk) + (1|sirepop),
                        ziformula =  ~  fst_wc_scaled + abs_tave_diff_midparent_sep_jul_scaled,
                        data = plants_gf_seeds,
                        family = nbinom2)
summary(seedsall.fst.noprecip.mod)
# temp and fst significant in zi part, fst marginal in conditional part

# does climate have same significance as in within-population plants when fst is not in model?
seedsall.fst.nofst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled +
                                abs_ppt_mm_diff_midparent_apr_jul_scaled
                                 + (1|site/blk) + (1|sirepop),
                                 ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
                                 data = plants_gf_seeds,
                                 family = nbinom2)
summary(seedsall.fst.nofst.mod)
# no, probably less ability to detect effects because midparent climate range is narrower in gene flow crosses. temp is significant in zi part only. 



# other lifestages: full models -------------------------------------------

# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

# generate scaled size covariate
plants_gf$nov_size_scaled = as.vector(scale(plants_gf$nov_size))
plants_gf$mar_size_scaled = as.vector(scale(plants_gf$mar_size))

germ.fst.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled + fst_wc_scaled 
                       + (1|site/blk) + (1|sirepop), 
                    family = binomial, 
                    data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fst.mod)
nrow(filter(plants_gf, !is.na(nov_germ)))
# 16236

novsize.fst.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled + fst_wc_scaled
                          + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fst.mod)
nrow(filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
# 4302

marsurv.fst.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + fst_wc_scaled + nov_size_scaled
                          + (1|site/blk) + (1|sirepop), 
                          family = binomial, 
                          data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fst.mod)
nrow(filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
# 4269

marsize.fst.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + fst_wc_scaled + nov_size_scaled
                          + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fst.mod)
nrow(filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
# 2591

fruitcount.fst.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled + mar_size_scaled
                             + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fst.mod)
nrow(filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
# 3055

seeds.fst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled + mar_size_scaled
                        + (1|site/blk) + (1|sirepop), 
                        family = nbinom2,
                        ziformula = ~ 1,
                        data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fst.mod)
nrow(filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
# 3055



# reduced predictors: fst only --------------------------------------------

germ.fstonly.mod = glmmTMB(nov_germ ~ fst_wc_scaled 
                           + (1|site/blk) + (1|sirepop), 
                       family = binomial, 
                       data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fstonly.mod)


novsize.fstonly.mod = glmmTMB(nov_size ~ fst_wc_scaled
                              + (1|site/blk) + (1|sirepop), 
                              family = gaussian,
                              data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fstonly.mod)


marsurv.fstonly.mod = glmmTMB(mar_surv ~ fst_wc_scaled + nov_size_scaled
                              + (1|site/blk) + (1|sirepop), 
                              family = binomial, 
                              data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fstonly.mod)


marsize.fstonly.mod = glmmTMB(mar_size ~ fst_wc_scaled + nov_size_scaled
                              + (1|site/blk) + (1|sirepop), 
                              family = gaussian,
                              data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fstonly.mod)


fruitcount.fstonly.mod = glmmTMB(fruit_count ~ fst_wc_scaled + mar_size_scaled
                                 + (1|site/blk) + (1|sirepop), 
                                 family = nbinom2,
                                 ziformula = ~ 1,
                                 data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fstonly.mod)


seeds.fstonly.mod = glmmTMB(total_est_seeds ~ fst_wc_scaled + mar_size_scaled
                            + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fstonly.mod)



# reduced predictors: temperature only ------------------------------------

germ.fsttemponly.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled 
                               + (1|site/blk) + (1|sirepop), 
                           family = binomial, 
                           data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fsttemponly.mod)


novsize.fsttemponly.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fsttemponly.mod)


marsurv.fsttemponly.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + nov_size_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = binomial, 
                          data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fsttemponly.mod)


marsize.fsttemponly.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + nov_size_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fsttemponly.mod)


fruitcount.fsttemponly.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + mar_size_scaled
                                     + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fsttemponly.mod)


seeds.fsttemponly.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + mar_size_scaled
                                + (1|site/blk) + (1|sirepop), 
                                family = nbinom2,
                                ziformula = ~ 1,
                                data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fsttemponly.mod)



# reduced predictors: precipitation only ----------------------------------

fruitcount.fstpptonly.mod = glmmTMB(fruit_count ~ abs_ppt_mm_diff_midparent_apr_jul_scaled + mar_size_scaled
                                    + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fstpptonly.mod)


seeds.fstpptonly.mod = glmmTMB(total_est_seeds ~ abs_ppt_mm_diff_midparent_apr_jul_scaled  + mar_size_scaled
                               + (1|site/blk) + (1|sirepop), 
                        family = nbinom2,
                        ziformula = ~ 1,
                        data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fstpptonly.mod)



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

results_fst = results_new %>% 
  unite(model, model, part, remove = TRUE, sep = ".") %>% 
  gather(var, value, c(3:11)) %>% 
  mutate(var = ifelse(var %in% c("abs_tave_diff_midparent_apr_jul_scaled", "abs_tave_diff_midparent_dec_mar_scaled", 
                                 "abs_tave_diff_midparent_sep_jul_scaled", "abs_tave_diff_midparent_sep_nov_scaled", 
                                 "abs_tave_diff_midparent_apr_may_scaled"), "temp", var)) %>%
  mutate(var = ifelse(var %in% c("abs_ppt_mm_diff_midparent_apr_jul_scaled", "abs_ppt_mm_diff_midparent_apr_may_scaled"), "ppt", var)) %>% 
  filter(!is.na(value)) %>%
  unite(variable, var, rowname, sep = ":") %>% 
  spread(variable, value)

results_round = rapply(object = results_fst, f = round, classes = "numeric", how = "replace", digits = 3) 

write.csv(results_round, "results/fst_table.csv", row.names = FALSE)




