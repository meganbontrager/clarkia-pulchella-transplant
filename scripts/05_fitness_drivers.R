# what drives fitness of within-pops groups?

library(tidyverse)
library(cowplot)
library(glmmTMB)
library(GGally)
library(ggeffects)

# coefficient plot code is from here:
# https://www.fromthebottomoftheheap.net/2017/05/04/compare-mgcv-with-glmmTMB/

# load the data ----------------------------------------------------

plants_wi = read.csv("data_processed/all_data.csv") %>% filter(type == "WI")
plants_wi$local_foreign = ifelse((plants_wi$sirepop == "AQ" & plants_wi$site == "AQ")|(plants_wi$sirepop == "AD" & plants_wi$site == "AD"), "local", "foreign")
plants = read.csv("data_processed/all_data.csv") 

# look at correlation of predictors --------------------------------

preds_temp = unique(plants_wi[,c("abs_tave_diff_sep_jul_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_sep_jul_scaled")])
ggpairs(preds_temp)
# ppt is correlated with temp and he, fis is corr with temp

preds_temp = unique(plants_wi[,c("abs_tave_diff_sep_nov_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_sep_nov_scaled")])
ggpairs(preds_temp)
# ppt is correlated with temp, fis

preds_temp = unique(plants_wi[,c("abs_tave_diff_dec_mar_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_dec_mar_scaled")])
ggpairs(preds_temp)
# ppt, fis are correlated with temp

preds_temp = unique(plants_wi[,c("abs_tave_diff_apr_may_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_apr_may_scaled")])
ggpairs(preds_temp)
# no correlations

preds_temp = unique(plants_wi[,c("abs_tave_diff_apr_jul_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_apr_jul_scaled")])
ggpairs(preds_temp)
# no correlations

preds_temp = unique(plants_wi[,c("abs_tave_diff_sep_jul_scaled", "fis_scaled", "he_scaled", "abs_ppt_mm_diff_apr_jul_scaled")])
ggpairs(preds_temp)
# no correlations


# total fitness, estimated by seeds --------------------------------

plants_wi_seeds = filter(plants_wi, !is.na(total_est_seeds))

# use zero inflated neagtive binomial to examine predictors of fitness
# with precipitation in spring/summer, when it is no longer correlated with temperature
seedsall.wi.mod = glmmTMB(total_est_seeds ~  abs_tave_diff_sep_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled
                          + (1|site/blk) + (1|sirepop), 
                                   ziformula = ~ abs_tave_diff_sep_jul_scaled,
                                   data = plants_wi_seeds, 
                                   family = nbinom2)
summary(seedsall.wi.mod) 


# other lifestages  -----------------------------------------------------
# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

plants_wi$nov_size_scaled = as.vector(scale(plants_wi$nov_size))
plants_wi$mar_size_scaled = as.vector(scale(plants_wi$mar_size))

germ.wi.mod = glmmTMB(nov_germ ~ abs_tave_diff_sep_nov_scaled 
                      + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = binomial, 
                      data = filter(plants_wi, !is.na(nov_germ)))
summary(germ.wi.mod)

novsize.wi.mod = glmmTMB(nov_size ~ abs_tave_diff_sep_nov_scaled
                         + (1|site/blk) + (1|sirepop), 
                         family = gaussian,
                         data = filter(plants_wi, nov_germ == 1, !is.na(nov_size)))
summary(novsize.wi.mod)

marsurv.wi.mod = glmmTMB(mar_surv ~ abs_tave_diff_dec_mar_scaled + nov_size_scaled
                         + (1|site) + (1|sirepop), 
                         family = binomial, 
                         data = filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.wi.mod)

marsize.wi.mod = glmmTMB(mar_size ~ abs_tave_diff_dec_mar_scaled  + nov_size_scaled 
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian,
                         data = filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.wi.mod)

fruitcount.wi.mod = glmmTMB(fruit_count ~ abs_tave_diff_apr_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled + mar_size_scaled 
                            + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.wi.mod)

seeds.wi.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_apr_jul_scaled  + abs_ppt_mm_diff_apr_jul_scaled + mar_size_scaled
                       + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.wi.mod)


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

results_wi = results_new %>% 
  unite(model, model, part, remove = TRUE, sep = ".") %>% 
  gather(var, value, c(3:10)) %>% 
  mutate(var = ifelse(var %in% c("abs_tave_diff_apr_jul_scaled", "abs_tave_diff_dec_mar_scaled", 
                                 "abs_tave_diff_sep_jul_scaled", "abs_tave_diff_sep_nov_scaled", 
                                 "abs_tave_diff_apr_may_scaled"), "temp", var)) %>%
  mutate(var = ifelse(var %in% c("abs_ppt_mm_diff_apr_jul_scaled", "abs_ppt_mm_diff_apr_may_scaled"), "ppt", var)) %>% 
  filter(!is.na(value)) %>%
  unite(variable, var, rowname, sep = ":") %>% 
  spread(variable, value)

results_round = rapply(object = results_wi, f = round, classes = "numeric", how = "replace", digits = 3) 

write.csv(results_round, "results/wi_table.csv", row.names = FALSE)


# plot of full model -----------------------------------------------------

mean(plants_wi_seeds$total_est_seeds)
pr.wi.seedsall.full = ggpredict(seedsall.wi.mod, terms = c("abs_tave_diff_sep_jul_scaled"), full.data = TRUE, typical = "mode")
plot(pr.wi.seedsall.full)
head(pr.wi.seedsall.full)
pr.wi.seedsall.full.sort = pr.wi.seedsall.full %>% arrange(x, observed)

distances = plants_wi_seeds %>% dplyr::select(x = abs_tave_diff_sep_jul_scaled, observed = total_est_seeds, dist, local_foreign, sire) %>% 
  arrange(x, observed)
pr.wi.seedsall.full.site = bind_cols(pr.wi.seedsall.full.sort, distances)

pr.wi.seedsall = ggpredict(seedsall.wi.mod, terms = c("abs_tave_diff_sep_jul_scaled"), typical = "mode")
plot(pr.wi.seedsall)

pr.wi.seedsall.full.means = pr.wi.seedsall.full.site %>% dplyr::group_by(x, local_foreign) %>% 
  dplyr::summarize(predicted = mean(predicted), conf.low = mean(conf.low), conf.high = mean(conf.high), dist = mean(dist))

plants_wi_means = plants_wi %>% group_by(abs_tave_diff_sep_jul_scaled, local_foreign) %>% summarize(x = mean(abs_tave_diff_sep_jul_scaled), total_est_seeds = mean(total_est_seeds, na.rm = TRUE))

sc = lm(plants$abs_tave_diff_sep_jul ~ 0 + plants$abs_tave_diff_sep_jul_scaled)
s = sc$coefficients[1]
m = mean(plants$abs_tave_diff_sep_jul)
names(plants_wi)


ggplot(data = pr.wi.seedsall, aes(x = x, y = predicted)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = pr.wi.seedsall.full.means, aes(x, predicted, shape = local_foreign, color = dist), size = 2.5) +
  geom_line(size = 1) + 
  # geom_point(data = plants_wi_means, aes(x, total_est_seeds), shape  = 18, size = 3) +
  scale_x_continuous(breaks = c((0-m)/s, (1-m)/s, (2-m)/s, (3-m)/s, (4-m)/s),
                     labels = c(0, 1, 2, 3, 4)) +
  xlab("Absolute temperature difference (°C)") +
  ylab("Predicted lifetime fitness") +
  guides(shape = FALSE, color = guide_colorbar(title = "Distance\nfrom gardens\n(km)", ticks = FALSE)) +
  theme(legend.title=element_text(size=12), legend.position = c(0.7, 0.7)) +
  scale_color_gradient(low = "dodgerblue", high  = "maroon2") +
  ylim(c(0, 47))

ggsave("figs/wi_temp.pdf", width = 4, height = 4)

# plot regression coefficients as panels ---------------------------
bTMB <- fixef(seedsall.wi.mod)$cond[-1]
seTMB <- diag(vcov(seedsall.wi.mod)$cond)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
df <- transform(df, upper = estimate + 1.96*sqrt(c(seTMB)), lower = estimate - 1.96*sqrt(c(seTMB)))

seeds.cond.reg = ggplot(df, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "",
       x = "Regression estimate on\nseed production given survival") +
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]))) +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(size = 10)) +
  xlim(c(-0.25, 0.25))

bTMB <- fixef(seedsall.wi.mod)$zi[-1]
seTMB <- diag(vcov(seedsall.wi.mod)$zi)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = (unname(c(bTMB))))
df <- transform(df, upper = (estimate + 1.96*sqrt(c(seTMB))), lower = (estimate - 1.96*sqrt(c(seTMB))))

seeds.zi.reg = ggplot(df,  aes(x = -estimate, y = term, xmax = -upper, xmin = -lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c(expression(T[diff]))) +
  labs(y = "",
       x = "Regression estimate on\nprobability of producing seeds") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(size = 10)) +
  xlim(c(-0.5, 0.5))

# plot whole thing

side.panels = plot_grid(seeds.cond.reg, seeds.zi.reg, ncol = 1, labels = c("B", "C"))

plot_grid(wi.main, side.panels, ncol = 2, labels = c("A", ""), rel_widths = c(1, 0.8))
ggsave("figs/wi_main.pdf", width = 7.5, height = 4.4)

# coefficient plots of within-pop performance ----------------------------------------

get_coef = function(mod){
  bTMB <- fixef(mod)$cond[-1]
  seTMB <- diag(vcov(mod)$cond)[-1]
  nms <- names(bTMB)
  df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
  df <- transform(df, upper = estimate + sqrt(c(seTMB)), lower = estimate - sqrt(c(seTMB)))
}

get_coef_zi = function(mod){
  bTMB <- fixef(mod)$zi[-1]
  seTMB <- diag(vcov(mod)$zi)[-1]
  nms <- names(bTMB)
  df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
  df <- transform(df, upper = estimate + sqrt(c(seTMB)), lower = estimate - sqrt(c(seTMB)))
}

germ.wi.coef = get_coef(germ.wi.mod)

germ.wi.reg = ggplot(germ.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\ngermination") +
  xlim(-0.31, 0.31) +
  theme(axis.title = element_text(size = 12))

novsize.wi.coef = get_coef(novsize.wi.mod)

novsize.wi.reg = ggplot(novsize.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nNovember size") +
  xlim(-0.2, 0.2) +
  theme(axis.title = element_text(size = 12))

marsurv.wi.coef = get_coef(marsurv.wi.mod)

marsurv.wi.reg = ggplot(marsurv.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]), "Nov.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nMarch survival") +
  xlim(-0.5, 0.5) +
  theme(axis.title = element_text(size = 12))

marsize.wi.coef = get_coef(marsize.wi.mod)

marsize.wi.reg = ggplot(marsize.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]), "Nov.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nMarch size") +
  xlim(-0.82, 0.82) +
  theme(axis.title = element_text(size = 12))

sumsurv.wi.coef = get_coef(sumsurv.wi.mod)

sumsurv.wi.reg = ggplot(sumsurv.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nsummer survival") +
  xlim(-0.53, 0.53) +
  theme(axis.title = element_text(size = 12))

fruitcount.wi.coef = get_coef(fruitcount.wi.mod)

fruitcount.wi.reg = ggplot(fruitcount.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nfruit production") +
  xlim(-0.57, 0.57) +
  theme(axis.title = element_text(size = 12))

seeds.wi.coef = get_coef(seeds.wi.mod)

seeds.wi.reg = ggplot(seeds.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +   
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nseed production") +
  xlim(-0.41, 0.41) +
  theme(axis.title = element_text(size = 12))

seedsall.wi.coef = get_coef(seedsall.wi.mod)

seedsall.wi.reg = ggplot(seedsall.wi.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nconditional seed production") +
  xlim(-0.2, 0.2) +
  theme(axis.title = element_text(size = 12))


seedsall.wi.coef.zi = get_coef_zi(seedsall.wi.mod)

seedsall.wi.reg.zi = ggplot(seedsall.wi.coef.zi, aes(x = -estimate, y = term, xmax = -upper, xmin = -lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nprobability of producing seeds") +
  xlim(-0.4, 0.4) +
  theme(axis.title = element_text(size = 12))


plot_grid(germ.wi.reg, novsize.wi.reg, marsurv.wi.reg,
          marsize.wi.reg, sumsurv.wi.reg, fruitcount.wi.reg,
          seeds.wi.reg, seedsall.wi.reg, seedsall.wi.reg.zi,
          ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
ggsave("figs/wi_estimates.pdf", height = 8, width = 8)




# plot survival curves -----------------
codes = read.csv("data_processed/pops_table.csv")

wi_curve = plants_wi %>% group_by(sirepop) %>% 
  summarise(start = n()/n(), nov_germ = mean(nov_germ), mar_surv = mean(mar_surv, na.rm = TRUE), fruits_any = mean(fruits_any, na.rm = TRUE)) %>% 
  gather("timepoint", "value", 2:5)

wi_curve = left_join(wi_curve, codes, by = c("sirepop" = "code"))

wi_curve$timepoint = ifelse(wi_curve$timepoint == "start", 1, 
                            ifelse(wi_curve$timepoint == "nov_germ", 2,
                                   ifelse(wi_curve$timepoint == "mar_surv", 3, 4)))

p1 = ggplot(wi_curve, aes(x = timepoint, y = value, color = Map.ID)) +
  geom_line()

gf_curve = plants_gf %>% group_by(sirepop) %>% 
  summarise(start = n()/n(), nov_germ = mean(nov_germ), mar_surv = mean(mar_surv, na.rm = TRUE), fruits_any = mean(fruits_any, na.rm = TRUE)) %>% 
  gather("timepoint", "value", 2:5)


gf_curve = left_join(gf_curve, codes, by = c("sirepop" = "code"))

gf_curve$timepoint = ifelse(gf_curve$timepoint == "start", 1, 
                            ifelse(gf_curve$timepoint == "nov_germ", 2,
                                   ifelse(gf_curve$timepoint == "mar_surv", 3, 4)))

p2 = ggplot(gf_curve, aes(x = timepoint, y = value, color = Map.ID)) +
  geom_line() +
  ylim(0,1)

plot_grid(p1, p2)

pops.wi.mod = glmmTMB(total_est_seeds ~  sirepop
                      + (1|site/blk) + (1|sire), 
                      ziformula = ~ sirepop,
                      data = plants_wi_seeds, 
                      family = nbinom2)
summary(pops.wi.mod) 

pred1 = ggpredict(pops.wi.mod, terms = "sirepop")
plot(pred1)

