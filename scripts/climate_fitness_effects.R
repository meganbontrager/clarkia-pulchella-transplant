# what drives fitness of within-pops groups?


# libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(glmmTMB) # version 0.2.2.0
packageVersion("glmmTMB")
library(GGally)
library(ggeffects) # version 0.3.0
packageVersion("ggeffects")



# load the data -----------------------------------------------------------

# load within-population crosses
plants_wi = read.csv("data/all_data.csv") %>% filter(type == "WI") %>% arrange(abs_tave_diff_sep_jul_scaled)
# 16125

# load full data frame (only used to rescale axis in plot)
plants = read.csv("data/all_data.csv") %>% arrange(abs_tave_diff_sep_jul_scaled)

# create local/foreign column (for plotting)
plants_wi$local_foreign = ifelse((plants_wi$sirepop == "AQ" & plants_wi$site == "AQ")|(plants_wi$sirepop == "AD" & plants_wi$site == "AD"), "local", "foreign")

# remove plants w/o reproduction estimates
plants_wi_seeds = filter(plants_wi, !is.na(total_est_seeds))
nrow(plants_wi_seeds)
# 13685



# look at correlation of predictors ---------------------------------------

preds_temp = unique(plants_wi[,c("abs_tave_diff_apr_jul_scaled", "abs_ppt_mm_diff_apr_jul_scaled")])
ggpairs(preds_temp)
# no strong correlation

preds_temp = unique(plants_wi[,c("abs_tave_diff_sep_jul_scaled", "abs_ppt_mm_diff_apr_jul_scaled")])
ggpairs(preds_temp)
# no strong correlation



# total fitness, estimated by seeds ---------------------------------------

# use zero inflated neagtive binomial to examine predictors of fitness
# with precipitation in spring/summer, when it is no longer correlated with temperature
seedsall.wi.mod = glmmTMB(total_est_seeds ~  abs_tave_diff_sep_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled
                          + (1|site/blk) + (1|sirepop), 
                          ziformula = ~ abs_tave_diff_sep_jul_scaled,
                          data = plants_wi_seeds, 
                          family = nbinom2)
summary(seedsall.wi.mod) 
nrow(plants_wi_seeds)
# 13685


# other lifestages  -------------------------------------------------------
# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

# scale size to use as a covariate
plants_wi$nov_size_scaled = as.vector(scale(plants_wi$nov_size))
plants_wi$mar_size_scaled = as.vector(scale(plants_wi$mar_size))

germ.wi.mod = glmmTMB(nov_germ ~ abs_tave_diff_sep_nov_scaled 
                      + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = binomial, 
                      data = filter(plants_wi, !is.na(nov_germ)))
summary(germ.wi.mod)
nrow(filter(plants_wi, !is.na(nov_germ)))
# 16125

novsize.wi.mod = glmmTMB(nov_size ~ abs_tave_diff_sep_nov_scaled
                         + (1|site/blk) + (1|sirepop), 
                         family = gaussian,
                         data = filter(plants_wi, nov_germ == 1, !is.na(nov_size)))
summary(novsize.wi.mod)
nrow(filter(plants_wi, nov_germ == 1, !is.na(nov_size)))
# 4303

marsurv.wi.mod = glmmTMB(mar_surv ~ abs_tave_diff_dec_mar_scaled + nov_size_scaled
                         + (1|site) + (1|sirepop), 
                         family = binomial, 
                         data = filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.wi.mod)
nrow(filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
# 4262

marsize.wi.mod = glmmTMB(mar_size ~ abs_tave_diff_dec_mar_scaled  + nov_size_scaled 
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian,
                         data = filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.wi.mod)
nrow(filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
# 2658

fruitcount.wi.mod = glmmTMB(fruit_count ~ abs_tave_diff_apr_jul_scaled + abs_ppt_mm_diff_apr_jul_scaled + mar_size_scaled 
                            + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.wi.mod)
nrow(filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
# 3176

seeds.wi.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_apr_jul_scaled  + abs_ppt_mm_diff_apr_jul_scaled + mar_size_scaled
                       + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.wi.mod)
nrow(filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
# 3176


# build a results table ---------------------------------------------------

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

# organize table
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

# round to 3 decimal pts
results_round = rapply(object = results_wi, f = round, classes = "numeric", how = "replace", digits = 3) 

# write.csv(results_round, "results/wi_table.csv", row.names = FALSE)


# plot of full model ------------------------------------------------------

# mean precip for plotting (should be close to 0 but not exactly since data has been subsetted since scaling)
mean(plants_wi_seeds$abs_ppt_mm_diff_apr_jul_scaled)
# -0.1238543

# mean of raw data
mean(plants_wi_seeds$total_est_seeds)
# 11.12115

# generate model predictions
pr.wi.seedsall = ggpredict(seedsall.wi.mod, terms = c("abs_tave_diff_sep_jul_scaled", "abs_ppt_mm_diff_apr_jul_scaled [-0.1238543]"), ci.lvl = 0.95, type = "fe", typical = "median", pretty = FALSE)
plot(pr.wi.seedsall)
mean(pr.wi.seedsall$predicted)

# generate raw means
plants_wi_means = plants_wi %>% group_by(abs_tave_diff_sep_jul_scaled, local_foreign, dist) %>% summarize(x = mean(abs_tave_diff_sep_jul_scaled), total_est_seeds = mean(total_est_seeds, na.rm = TRUE))

# rescale x-axis to original units
sc = lm(plants$abs_tave_diff_sep_jul ~ 0 + plants$abs_tave_diff_sep_jul_scaled)
s = sc$coefficients[1]
m = mean(plants$abs_tave_diff_sep_jul)

# plot
ggplot(data = pr.wi.seedsall, aes(x = x, y = predicted)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(size = 1) + 
  geom_point(data = plants_wi_means, aes(x, total_est_seeds, shape = local_foreign, color = dist), size = 3) +
  scale_x_continuous(breaks = c((0-m)/s, (1-m)/s, (2-m)/s, (3-m)/s, (4-m)/s),
                     labels = c(0, 1, 2, 3, 4)) +
  xlab("Absolute temperature difference (°C)") +
  ylab("Predicted lifetime fitness") +
  guides(shape = FALSE, color = guide_colorbar(title = "Distance\nfrom gardens\n(km)", ticks = FALSE)) +
  theme(legend.title=element_text(size=12), legend.position = c(0.7, 0.7)) +
  scale_color_gradient(low = "dodgerblue", high  = "maroon2") +
  ylim(c(0, 40))


# plot regression coefficients as panels ----------------------------------

# coefficient plot code is from here:
# https://www.fromthebottomoftheheap.net/2017/05/04/compare-mgcv-with-glmmTMB/

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

# plot the whole thing
plot_grid(seeds.cond.reg, seeds.zi.reg, ncol = 1)
