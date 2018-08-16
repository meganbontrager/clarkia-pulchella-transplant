# are there effects of gene flow in addition to environment?

library(tidyverse)
library(cowplot)
library(glmmTMB)
library(GGally)
library(ggeffects)
library(lemon)

# load data ---------------------------------------------------------------

plants = read.csv("data_processed/all_data.csv") 
plants_wi = read.csv("data_processed/all_data.csv") %>% filter(type == "WI")
plants_gf = read.csv("data_processed/all_data.csv") %>% filter(type == "GF")



# do ranges look as expected? ---------------------------------------------
hist(plants_wi$abs_tave_diff_midparent_sep_jul)
hist(plants_gf$abs_tave_diff_midparent_sep_jul)
hist(plants_wi$abs_tave_diff_midparent_sep_jul_scaled)
hist(plants_gf$abs_tave_diff_midparent_sep_jul_scaled)

plants$local_foreign = ifelse((plants$sirepop == "AQ" & plants$site == "AQ")|(plants$sirepop == "AD" & plants$site == "AD"), "local", "foreign")

plants_for_hist = plants %>% select(type, local_foreign, abs_ppt_mm_diff_midparent_apr_jul, abs_tave_diff_midparent_sep_jul, sirepop, dampop) %>% distinct()

p1 = ggplot(plants_for_hist, aes(abs_ppt_mm_diff_midparent_apr_jul, fill = local_foreign)) +
  geom_histogram(bins = 12) +
  scale_fill_manual(values = c("grey40", "black")) +
  facet_rep_grid(type~.) +
  guides(fill = FALSE) +
  ylab("Number of crossing groups") +
  xlab("Absolute midparent\nprecipitation difference (mm)") +
  theme(strip.background = element_blank(), strip.text = element_blank(), axis.line=element_line())

p2 = ggplot(plants_for_hist, aes(abs_tave_diff_midparent_sep_jul, fill = local_foreign)) +
  geom_histogram(bins = 12) +
  scale_fill_manual(values = c("grey40", "black")) +
  facet_rep_grid(type~.) +
  guides(fill = FALSE) +
  ylab("Number of crossing groups") +
  xlab("Absolute midparent\ntemperature difference (°C)") +
  theme(strip.background = element_blank(), strip.text = element_blank(), axis.line=element_line())

c= plot_grid(p1, p2, labels = c("A", "B"))
  
ggdraw(c) +
  draw_label("Between-\npopulations", 0.42, 0.95) +
  draw_label("Between-\npopulations", 0.92, 0.95) + 
  draw_label("Within-\npopulations", 0.42, 0.50) + 
  draw_label("Within-\npopulations", 0.92, 0.50)


# overall fitness ---------------------------------------------------------

plants_seeds = plants %>% filter(!is.na(total_est_seeds))

seedsall.gf.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled + type + abs_ppt_mm_diff_midparent_apr_jul_scaled
                          + (1|site/blk) + (1|sirepop), 
                          ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled + type,
                          data = plants_seeds, 
                          family = nbinom2)
summary(seedsall.gf.mod) 

seedsall.gf.mod.int = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled * type
                          + (1|site/blk) + (1|sirepop), 
                          ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled * type,
                          data = plants_seeds, 
                          family = nbinom2)
summary(seedsall.gf.mod.int) 


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


AIC(seedsall.gf.mod, seedsall.gf.noppt.mod, seedsall.gf.notype.mod)

seedsall.gf.simplezi.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_sep_jul_scaled +
                                   abs_ppt_mm_diff_midparent_apr_jul_scaled + type
                                   + (1|site/blk) + (1|sirepop),
                                   ziformula =  ~ abs_tave_diff_midparent_sep_jul_scaled,
                                   data = plants_seeds,
                                   family = nbinom2)
summary(seedsall.gf.simplezi.mod)


# other lifestages  -----------------------------------------------------
# using size as a covariate, seasonal temperature, and datasets that are filtered to only plants alive during the last census

plants$nov_size_scaled = as.vector(scale(plants$nov_size))
plants$mar_size_scaled = as.vector(scale(plants$mar_size))

germ.gf.mod = glmmTMB(nov_germ ~ abs_tave_diff_midparent_sep_nov_scaled + type
                      + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = binomial, 
                      data = filter(plants, !is.na(nov_germ)))
summary(germ.gf.mod)

novsize.gf.mod = glmmTMB(nov_size ~ abs_tave_diff_midparent_sep_nov_scaled + type
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = gaussian,
                      data = filter(plants, nov_germ == 1, !is.na(nov_size)))
summary(novsize.gf.mod)

marsurv.gf.mod = glmmTMB(mar_surv ~ abs_tave_diff_midparent_dec_mar_scaled + type + nov_size_scaled
                         + (1|site) + (1|sirepop), 
                       family = binomial, 
                       data = filter(plants, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.gf.mod)

marsize.gf.mod = glmmTMB(mar_size ~ abs_tave_diff_midparent_dec_mar_scaled + type + nov_size_scaled
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = gaussian,
                      data = filter(plants, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.gf.mod)

fruitcount.gf.mod = glmmTMB(fruit_count ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + type + mar_size_scaled
                            + (1|site/blk) + (1|sirepop), 
                            family = nbinom2,
                            ziformula = ~ 1,
                            data = filter(plants, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.gf.mod)

seeds.gf.mod = glmmTMB(total_est_seeds ~ abs_tave_diff_midparent_apr_jul_scaled + abs_ppt_mm_diff_midparent_apr_jul_scaled + type + mar_size_scaled
                       + (1|site/blk) + (1|sirepop), 
                       family = nbinom2,
                       ziformula = ~ 1,
                       data = filter(plants, sum_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.gf.mod)


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

write.csv(results_round, "results/gf_table.csv", row.names = FALSE)


# plot of full model -----------------------------------------------------

pr.gf.seedsall = ggpredict(seedsall.gf.mod, c("abs_tave_diff_midparent_sep_jul_scaled", "type"))
plot(pr.gf.seedsall)

pr.gf.seedsall2 = ggpredict(seedsall.gf.mod, c("abs_ppt_mm_diff_midparent_apr_jul_scaled", "type"))
plot(pr.gf.seedsall2)

pr.gf.seedsall3 = ggpredict(seedsall.gf.noppt.mod, c("type"))
plot(pr.gf.seedsall3)


gf_max = max(plants_gf$abs_tave_diff_midparent_sep_jul_scaled)
gf_min = min(plants_gf$abs_tave_diff_midparent_sep_jul_scaled)
pr.gf.seedsall.wi = pr.gf.seedsall %>% filter(group == "WI")
pr.gf.seedsall.gf = pr.gf.seedsall %>% filter(group == "GF") %>% filter(x <= gf_max, x >= gf_min)

pr.gf.seedsall.noppt = ggpredict(seedsall.gf.noppt.mod, c("abs_tave_diff_midparent_sep_jul_scaled", "type"))
plot(pr.gf.seedsall.noppt)

pr.gf.frcount = ggpredict(fruitcount.gf.mod, c("abs_tave_diff_midparent_apr_jul_scaled", "type"))
plot(pr.gf.frcount)

# plants_gf_means = plants_seeds %>% dplyr::group_by(abs_tave_diff_midparent_sep_jul_scaled, type) %>% 
#   dplyr::summarize(total_est_seeds = mean(total_est_seeds))

sc = lm(plants$abs_tave_diff_midparent_sep_jul ~ 0 + plants$abs_tave_diff_midparent_sep_jul_scaled)
s = sc$coefficients[1]
m = mean(plants$abs_tave_diff_midparent_sep_jul)

# gf.main = 
  ggplot() +
  geom_line(data = pr.gf.seedsall.wi, aes(x, predicted), color = "grey40", size = 1) + 
  geom_ribbon(data = pr.gf.seedsall.wi, aes(x= x, ymin = conf.low, ymax = conf.high), fill = "grey40", alpha = 0.3) +
  geom_line(data = pr.gf.seedsall.gf, aes(x, predicted), color = "green4", size = 1) + 
  geom_ribbon(data = pr.gf.seedsall.gf, aes(x= x, ymin = conf.low, ymax = conf.high), fill = "green4", alpha = 0.3) +
  # geom_point(data = plants_gf_means, aes(x = abs_tave_diff_midparent_sep_jul_scaled, y = total_est_seeds, color = factor(type)), alpha = 0.4) +
  scale_x_continuous(breaks = c((0-m)/s, (1-m)/s, (2-m)/s, (3-m)/s, (4-m)/s), 
                     labels = c(0, 1, 2, 3, 4)) +
  xlab("Absolute midparent difference from\nhistoric average temperature (°C)") +
  ylab("Predicted lifetime fitness") +
  ylim(c(0,51))

pr.gf.seeds = ggpredict(seeds.gf.mod, c("abs_tave_diff_midparent_apr_jul_scaled", "type"))
pr.gf.seeds.wi = pr.gf.seeds %>% filter(group == "WI")
pr.gf.seeds.gf = pr.gf.seeds %>% filter(group == "GF") %>% filter(x <= gf_max, x >= gf_min)

# plants_gf_means2 = plants_seeds %>% filter(mar_surv == 1) %>% dplyr::group_by(abs_tave_diff_midparent_apr_jul_scaled, type) %>% 
#   dplyr::summarize(total_est_seeds = mean(total_est_seeds))

# gf.main2 = 
  ggplot() +
  geom_line(data = pr.gf.seeds.wi, aes(x, predicted), color = "grey40", size = 1) + 
  geom_ribbon(data = pr.gf.seeds.wi, aes(x= x, ymin = conf.low, ymax = conf.high), fill = "grey40", alpha = 0.3) +
  geom_line(data = pr.gf.seeds.gf, aes(x, predicted), color = "green3", size = 1) + 
  geom_ribbon(data = pr.gf.seeds.gf, aes(x= x, ymin = conf.low, ymax = conf.high), fill = "green3", alpha = 0.3) +
  # geom_point(data = plants_gf_means2, aes(x = abs_tave_diff_midparent_apr_jul_scaled, y = total_est_seeds, color = factor(type)), alpha = 0.4) +
  scale_color_manual(values = c("red", "blue")) +
  scale_x_continuous(breaks = c((0-m)/s, (1-m)/s, (2-m)/s, (3-m)/s, (4-m)/s), 
                     labels = c(0, 1, 2, 3, 4)) +
  xlab("Absolute midparent\ntemperature difference (°C)") +
  ylab("Predicted per capita\nseed production") +
  ylim(0,250) +
  guides(color = FALSE)

# plot regression coefficients as panels
bTMB <- fixef(seedsall.gf.mod)$cond[-1]
seTMB <- diag(vcov(seedsall.gf.mod)$cond)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = unname(c(bTMB)))
df <- transform(df, upper = estimate + sqrt(c(seTMB)), lower = estimate - sqrt(c(seTMB)))
df[2,c(2:4)] = -df[2,c(2:4)] # change wi to gf

seeds.cond.re = ggplot(df, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "",
       x = "Regression estimate on\nseed production given survival") +
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), "GF")) +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(size = 10)) +
  xlim(c(-0.2, 0.2))

bTMB <- fixef(seedsall.gf.mod)$zi[-1]
seTMB <- diag(vcov(seedsall.gf.mod)$zi)[-1]
nms <- names(bTMB)
df <- data.frame(term = c(nms), estimate = (unname(c(bTMB))))
df <- transform(df, upper = (estimate + sqrt(c(seTMB))), lower = (estimate - sqrt(c(seTMB))))
df[2,c(2:4)] = -df[2,c(2:4)]

seeds.zi.reg = ggplot(df,  aes(x = -estimate, y = term, xmax = -upper, xmin = -lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c(expression(T[diff]), "GF")) +
  labs(y = "",
       x = "Regression estimate on\nprobability of producing seeds") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(size = 10)) +
  xlim(c(-0.3, 0.3))

# plot whole thing

side.panels = plot_grid(seeds.cond.reg, seeds.zi.reg, ncol = 1, labels = c("B", "C"))

plot_grid(gf.main, side.panels, ncol = 2, labels = c("A", ""), rel_gfdths = c(1, 0.8))
ggsave("figs/gf_main.pdf", width = 7.5, height = 4.4)

# coefficient plots of gene flow and performance ----------------------------------------
# needs updating!

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

germ.gf.coef = get_coef(germ.gf.mod)

germ.gf.reg = ggplot(germ.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\ngermination") +
  xlim(-0.31, 0.31) +
  theme(axis.title = element_text(size = 12))

novsize.gf.coef = get_coef(novsize.gf.mod)

novsize.gf.reg = ggplot(novsize.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nNovember size") +
  xlim(-0.2, 0.2) +
  theme(axis.title = element_text(size = 12))

marsurv.gf.coef = get_coef(marsurv.gf.mod)

marsurv.gf.reg = ggplot(marsurv.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]), "Nov.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nMarch survival") +
  xlim(-0.5, 0.5) +
  theme(axis.title = element_text(size = 12))

marsize.gf.coef = get_coef(marsize.gf.mod)

marsize.gf.reg = ggplot(marsize.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]), "Nov.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nMarch size") +
  xlim(-0.82, 0.82) +
  theme(axis.title = element_text(size = 12))

sumsurv.gf.coef = get_coef(sumsurv.gf.mod)

sumsurv.gf.reg = ggplot(sumsurv.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nsummer survival") +
  xlim(-0.53, 0.53) +
  theme(axis.title = element_text(size = 12))

fruitcount.gf.coef = get_coef(fruitcount.gf.mod)

fruitcount.gf.reg = ggplot(fruitcount.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nfruit production") +
  xlim(-0.57, 0.57) +
  theme(axis.title = element_text(size = 12))

seeds.gf.coef = get_coef(seeds.gf.mod)

seeds.gf.reg = ggplot(seeds.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +   
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]), "Mar.\nsize")) +
  labs(y = "Parameter",
       x = "Regression estimate on\nseed production") +
  xlim(-0.41, 0.41) +
  theme(axis.title = element_text(size = 12))

seedsall.gf.coef = get_coef(seedsall.gf.mod)

seedsall.gf.reg = ggplot(seedsall.gf.coef, aes(x = estimate, y = term, xmax = upper, xmin = lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(P[diff]), expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nconditional seed production") +
  xlim(-0.2, 0.2) +
  theme(axis.title = element_text(size = 12))


seedsall.gf.coef.zi = get_coef_zi(seedsall.gf.mod)

seedsall.gf.reg.zi = ggplot(seedsall.gf.coef.zi, aes(x = -estimate, y = term, xmax = -upper, xmin = -lower)) +
  geom_point() +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_discrete(labels = c(expression(T[diff]), expression(F[IS]), expression(H[E]))) +
  labs(y = "Parameter",
       x = "Regression estimate on\nprobability of producing seeds") +
  xlim(-0.4, 0.4) +
  theme(axis.title = element_text(size = 12))


plot_grid(germ.gf.reg, novsize.gf.reg, marsurv.gf.reg,
          marsize.gf.reg, sumsurv.gf.reg, fruitcount.gf.reg,
          seeds.gf.reg, seedsall.gf.reg, seedsall.gf.reg.zi,
          ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
ggsave("figs/wi_estimates.pdf", height = 8, width = 8)


