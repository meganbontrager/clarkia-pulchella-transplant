library(tidyverse)
library(cowplot)
library(glmmTMB)
library(GGally)
library(ggeffects)

plants_gf = read.csv("data_processed/all_data.csv") %>% filter(type == "GF")

# overall fitness ---------------------------------------------------------

# filter to rows with seeds
plants_gf_seeds = plants_gf %>% filter(!is.na(total_est_seeds))

# correlation of predictors -----------------------------------------------

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_sep_jul_scaled", "fst_wc_scaled", "abs_ppt_mm_diff_midparent_apr_jul_scaled")])
ggpairs(preds_temp)
# ppt is correlated with fst at 0.635

preds_temp = unique(plants_gf[,c("ppt_mm_diff_midparent_apr_jul_scaled", "fst_wc_scaled")])
ggpairs(preds_temp) # 0.603

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_sep_nov_scaled", "fst_wc_scaled")])
ggpairs(preds_temp)

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_dec_mar_scaled", "fst_wc_scaled")])
ggpairs(preds_temp)

preds_temp = unique(plants_gf[,c("abs_tave_diff_midparent_apr_jul_scaled", "fst_wc_scaled", "abs_ppt_mm_diff_midparent_apr_jul_scaled")])
ggpairs(preds_temp)

seedsall.fst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + 
                          abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled
                                    + (1|site/blk) + (1|sirepop), 
                                    ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + fst_wc_scaled,
                                    data = plants_gf_seeds, 
                                    family = nbinom2)
summary(seedsall.fst.mod) 

# precip is not significant in univariate models (below) and is co-linear with fst (high fst == high difference in precip), so trying without it.

# seedsall.fst.noprecip.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + fst_wc_scaled
#                                  + (1|site/blk) + (1|sirepop), 
#                         ziformula =  ~  fst_wc_scaled + abs_tave_diff_midparent_sep_jul_scaled,
#                         data = plants_gf_seeds, 
#                         family = nbinom2)
# summary(seedsall.fst.noprecip.mod) 
 
# above won't converge. univariate models indicate that fst and temp matter for zi, but only fst for cond. trying that:
# seedsall.fst.noprecip2.mod = glmmTMB(total_est_seeds ~ fst_wc_scaled
#                                   + (1|site/blk) + (1|sirepop), 
#                                  ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + fst_wc_scaled,
#                                  data = plants_gf_seeds, 
#                                  family = nbinom2)
# summary(seedsall.fst.noprecip2.mod) 
 

# seedsall.fst.nofst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled +
#                                 abs_ppt_mm_diff_midparent_apr_jul_scaled
#                                  + (1|site/blk) + (1|sirepop),
#                                  ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
#                                  data = plants_gf_seeds,
#                                  family = nbinom2)
# summary(seedsall.fst.nofst.mod)

# univariate models

seedsall.fst.temponly.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled
                                    + (1|site/blk) + (1|sirepop), 
                                 ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
                                 data = plants_gf_seeds, 
                                 family = nbinom2)
summary(seedsall.fst.temponly.mod) 
# t sig
# conv w big re

seedsall.fst.fstonly.mod = glmmTMB(total_est_seeds ~ fst_wc_scaled
                                   + (1|site/blk) + (1|sirepop), 
                                 ziformula =  ~ fst_wc_scaled ,
                                 data = plants_gf_seeds, 
                                 family = nbinom2)
summary(seedsall.fst.fstonly.mod) 
# fst sig
# no conv w big re

seedsall.fst.preciponly.mod = glmmTMB(total_est_seeds ~ abs_ppt_mm_diff_midparent_apr_jul_scaled
                                      + (1|site/blk) + (1|sirepop), 
                                 ziformula =  ~ 1,
                                 data = plants_gf_seeds, 
                                 family = nbinom2)
summary(seedsall.fst.preciponly.mod) 
# ppt marg
# no conv w big re


# other lifestages  -----------------------------------------------------
# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

plants_gf$nov_size_scaled = as.vector(scale(plants_gf$nov_size))
plants_gf$mar_size_scaled = as.vector(scale(plants_gf$mar_size))

germ.fst.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled + fst_wc_scaled 
                       + (1|site/blk) + (1|sirepop), 
                    family = binomial, 
                    data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fst.mod)
# fst sig
# no conv w big re

germ.fstonly.mod = glmmTMB(nov_germ ~ fst_wc_scaled 
                           + (1|site/blk) + (1|sirepop), 
                       family = binomial, 
                       data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fstonly.mod)
# fst sig
# conv w big re

germ.fsttemponly.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled 
                               + (1|site/blk) + (1|sirepop), 
                           family = binomial, 
                           data = filter(plants_gf, !is.na(nov_germ)))
summary(germ.fsttemponly.mod)
# t ns
# conv w big re

novsize.fst.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled + fst_wc_scaled
                          + (1|site/blk) + (1|sirepop), 
                      family = gaussian,
                      data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fst.mod)
# ns
# no conv w big re

novsize.fstonly.mod = glmmTMB(nov_size ~ fst_wc_scaled
                              + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fstonly.mod)
# ns
# conv w big re

novsize.fsttemponly.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, nov_germ == 1, !is.na(nov_size)))
summary(novsize.fsttemponly.mod)
# marg
# conv w big re

marsurv.fst.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + fst_wc_scaled + nov_size_scaled
                          + (1|site/blk) + (1|sirepop), 
                       family = binomial, 
                       data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fst.mod)
# ns
# no conv w big re

marsurv.fstonly.mod = glmmTMB(mar_surv ~ fst_wc_scaled + nov_size_scaled
                              + (1|site/blk) + (1|sirepop), 
                          family = binomial, 
                          data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fstonly.mod)
# ns
# no conv w big re

marsurv.fsttemponly.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + nov_size_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = binomial, 
                          data = filter(plants_gf, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.fsttemponly.mod)
# ns
# no conv w big re


marsize.fst.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + fst_wc_scaled + nov_size_scaled
                          + (1|site/blk) + (1|sirepop), 
                      family = gaussian,
                      data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fst.mod)
# ns
# conv w big re

marsize.fstonly.mod = glmmTMB(mar_size ~ fst_wc_scaled + nov_size_scaled
                              + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fstonly.mod)
# sig
# no conv w big re

marsize.fsttemponly.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + nov_size_scaled
                                  + (1|site/blk) + (1|sirepop), 
                          family = gaussian,
                          data = filter(plants_gf, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.fsttemponly.mod)
# sig
# conv w big re

fruitcount.fst.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled + mar_size_scaled
                             + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fst.mod)
# ns
# conv w big re

fruitcount.fstonly.mod = glmmTMB(fruit_count ~ fst_wc_scaled + mar_size_scaled
                                 + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fstonly.mod)
# sig
# no conv w big re

fruitcount.fsttemponly.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + mar_size_scaled
                                     + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fsttemponly.mod)
# marg
# no conv w big re

fruitcount.fstpptonly.mod = glmmTMB(fruit_count ~ abs_ppt_mm_diff_midparent_apr_jul_scaled + mar_size_scaled
                                    + (1|site/blk) + (1|sirepop), 
                             family = nbinom2,
                             ziformula = ~ 1,
                             data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.fstpptonly.mod)
# ns
# conv w big re

seeds.fst.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + fst_wc_scaled + mar_size_scaled
                        + (1|site/blk) + (1|sirepop), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fst.mod)
# fst marg
# no conv w big re

seeds.fsttemponly.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + mar_size_scaled
                                + (1|site/blk) + (1|sirepop), 
                        family = nbinom2,
                        ziformula = ~ 1,
                        data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fsttemponly.mod)
# ns
# conv w big re

seeds.fstpptonly.mod = glmmTMB(total_est_seeds ~ abs_ppt_mm_diff_midparent_apr_jul_scaled  + mar_size_scaled
                               + (1|site/blk) + (1|sirepop), 
                        family = nbinom2,
                        ziformula = ~ 1,
                        data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fstpptonly.mod)
# ns
# no conv w big re

seeds.fstonly.mod = glmmTMB(total_est_seeds ~ fst_wc_scaled + mar_size_scaled
                            + (1|site/blk) + (1|sirepop), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants_gf, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.fstonly.mod)
# ns
# no conv w big re


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
  gather(var, value, c(3:13)) %>% 
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


pr.fst.seedsall.full = ggpredict(seedsall.fst.mod, c("fst_wc_scaled"), full.data = TRUE, typical = "mode")
plot(pr.fst.seedsall.full)

summary(pr.fst.seedsall.full$x)

pr.fst.seedsall.full.sort = pr.fst.seedsall.full %>% arrange(x, observed)

distances = plants_gf_seeds %>% dplyr::select(x = fst_wc_scaled, observed = total_est_seeds, dist, sire) %>% 
  arrange(x, observed)
summary(distances$x)

pr.fst.seedsall.full.site = bind_cols(pr.fst.seedsall.full.sort, distances)

plants_gf_means = plants_gf_seeds %>% dplyr::group_by(fst_wc_scaled) %>% 
  dplyr::summarize(total_est_seeds = mean(total_est_seeds))

pr.fst.seedsall = ggpredict(seedsall.fst.mod, c("fst_wc_scaled"), typical = "mode")
plot(pr.fst.seedsall)

pr.fst.seedsall.full.means = pr.fst.seedsall.full.site %>% dplyr::group_by(x) %>% 
  dplyr::summarize(predicted = mean(predicted), conf.low = mean(conf.low), conf.high = mean(conf.high), dist = mean(dist))

# plants_wi_means_local = plants_wi_seeds %>% dplyr::filter((sirepop == "AQ" & site == "AQ")|(sirepop == "AD" & site == "AD")) %>% dplyr::group_by(abs_tave_diff_sep_jul_scaled) %>% 
#   dplyr::summarize(total_est_seeds = mean(total_est_seeds)) 

sc = lm(plants$fst_wc_hfs ~ 0 + plants$fst_wc_scaled)
s = sc$coefficients[1]
m = mean(plants$fst_wc_hfs, na.rm = TRUE)

ggplot(data = pr.fst.seedsall, aes(x = x, y = predicted)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = pr.fst.seedsall.full.means, aes(x, predicted, color = dist), size = 2.5) +
  geom_line(size = 1) + 
  # geom_point(data = plants_wi_means, aes(x, total_est_seeds), shape  = 18, size = 3) +
  scale_x_continuous(breaks = c((0-m)/s, (0.05-m)/s, (0.1-m)/s, (0.15-m)/s, (0.2-m)/s, (0.25-m)/s),
                     labels = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
  xlab(expression(Genetic~distance~(F[ST]))) +
  ylab("Predicted lifetime fitness") +
  guides(shape = FALSE, color = guide_colorbar(title = "Distance\nfrom gardens\n(km)", ticks = FALSE)) +
  theme(legend.title=element_text(size=12), legend.position = c(0.05, 0.75)) +
  scale_color_gradient(low = "dodgerblue", high  = "maroon2") +
  ylim(c(0, 52))

ggsave("figs/fst_temp.pdf", width = 4, height = 4)

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
  xlim(c(-0.15, 0.15))
