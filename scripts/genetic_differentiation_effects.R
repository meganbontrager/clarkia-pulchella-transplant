library(tidyverse)
library(cowplot)
library(GGally)
library(glmmTMB) # version 0.2.2.0
library(ggeffects) # version 0.3.0 


# load data ---------------------------------------------------------------

# load and filter to between-population crosses
plants_gf = read.csv("data/all_data.csv") %>% filter(type == "GF")
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





# plot --------------------------------------------------------------------

# pr.fst.seedsall.full.site = bind_cols(pr.fst.seedsall.full.sort, distances)
# generate raw means
plants_gf_means = plants_gf_seeds %>% dplyr::group_by(fst_wc_scaled, fst_wc_hfs, dist, abs_tave_diff_sep_jul_scaled) %>% 
  dplyr::summarize(total_est_seeds = mean(total_est_seeds))

# generate model prediction and CI
pr.fst.seedsall = ggpredict(seedsall.fst.mod, c("fst_wc_scaled"), typical = "median")
mean(pr.fst.seedsall$predicted)
plot(pr.fst.seedsall)

# rescale x-axis back to original units
sc = lm(plants_gf$fst_wc_hfs ~ 0 + plants_gf$fst_wc_scaled)
s = sc$coefficients[1]
m = mean(plants_gf$fst_wc_hfs, na.rm = TRUE)

# main plot

ggplot(data = pr.fst.seedsall, aes(x = x, y = predicted)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(size = 1) + 
  geom_point(data = plants_gf_means, aes(fst_wc_scaled, total_est_seeds, color = dist), size = 3) +
  scale_x_continuous(breaks = c((0-m)/s, (0.05-m)/s, (0.1-m)/s, (0.15-m)/s, (0.2-m)/s, (0.25-m)/s),
                     labels = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
  xlab(expression(Genetic~differentiation~(F[ST]))) +
  ylab("Predicted lifetime fitness") +
  guides(shape = FALSE, color = guide_colorbar(title = "Distance\nfrom gardens\n(km)", ticks = FALSE)) +
  theme(legend.title=element_text(size=12), legend.position = c(0.05, 0.75)) +
  scale_color_gradient(low = "dodgerblue", high  = "maroon2") +
  ylim(c(0, 30))

# plot regression coefficients as panels

bTMB <- fixef(seedsall.fst.mod)$cond[-1]
seTMB <- diag(vcov(seedsall.fst.mod)$cond)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
df <- transform(df, upper = estimate + sqrt(c(seTMB)), lower = estimate - sqrt(c(seTMB)))

fst.cond.re = ggplot(df, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "",
       x = "Regression estimate on\nseed production given survival") +
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[ST]))) +
  theme(axis.title = element_text(size = 12)) +
  xlim(c(-0.15, 0.15))

bTMB <- fixef(seedsall.fst.mod)$zi[-1]
seTMB <- diag(vcov(seedsall.fst.mod)$zi)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
df <- transform(df, upper = estimate + sqrt(c(seTMB)), lower = estimate - sqrt(c(seTMB)))

fst.zi.re = ggplot(df, aes(x = -estimate, y = term, xmax = -upper, xmin = -lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "",
       x = "Regression estimate on\nprobability of producing seeds") +
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[ST]))) +
  theme(axis.title = element_text(size = 12)) +
  xlim(c(-0.17, 0.17))

plot_grid(fst.zi.re, fst.cond.re, ncol = 1, rel_heights = c(0.4, 0.6))


