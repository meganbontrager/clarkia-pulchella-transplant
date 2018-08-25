# do local edge populations outperform foreign populations?

# libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(glmmTMB) # version 0.2.2.0
packageVersion("glmmTMB")
library(ggeffects) # version 0.3.3
packageVersion("ggeffects")



# load data ---------------------------------------------------------------

# load data and filter to within-population crosses
plants_wi =  read.csv("data/all_data.csv") %>% filter(type == "WI")

# create local/foreign column where focal pops in their home sites are designated as local
plants_wi$local_foreign = ifelse((plants_wi$sirepop == "AQ" & plants_wi$site == "AQ")|(plants_wi$sirepop == "AD" & plants_wi$site == "AD"), "local", "foreign")

# filter to rows that have total estimated seeds (i.e. omit plants that were killed by gophers, etc.)
plants_wi_seeds = filter(plants_wi, !is.na(total_est_seeds))



# total fitness, estimated by seeds ---------------------------------------

# use zero inflated negtive binomial to compare average performance of locals v. foreigns
seedsall.lf.mod = glmmTMB(total_est_seeds ~ local_foreign 
                          + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                          ziformula = ~ local_foreign,
                          data = plants_wi_seeds, 
                          family = nbinom2)
summary(seedsall.lf.mod)
# no significant difference between locals and foreigns



# other lifestages with size as covariate, truncated datasets -------------

# scale size to use as a covariate
plants_wi$nov_size_scaled = as.vector(scale(plants_wi$nov_size))
plants_wi$mar_size_scaled = as.vector(scale(plants_wi$mar_size))

# germination
germ.lf.mod = glmmTMB(nov_germ ~ local_foreign 
                      + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                      family = binomial, 
                      data = plants_wi)
summary(germ.lf.mod)
nrow(plants_wi)
# 16125

# size after germination
novsize.lf.mod = glmmTMB(nov_size ~ local_foreign
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian, 
                         data = filter(plants_wi, !is.na(nov_size_scaled)))
summary(novsize.lf.mod)
nrow(filter(plants_wi, !is.na(nov_size_scaled)))
# 4303

# overwinter survival
marsurv.lf.mod = glmmTMB(mar_surv ~ local_foreign + nov_size_scaled
                         + (1|sirepop) + (1|site), 
                         family = binomial, 
                         data = filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
summary(marsurv.lf.mod)
nrow(filter(plants_wi, nov_germ == 1, !is.na(mar_surv), !is.na(nov_size_scaled)))
# 4262

# march size
marsize.lf.mod = glmmTMB(mar_size ~ local_foreign + nov_size_scaled
                         + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                         family = gaussian,
                         data = filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
summary(marsize.lf.mod)
nrow(filter(plants_wi, mar_surv == 1, !is.na(mar_size), !is.na(nov_size_scaled)))
# 2658

# fruit number
fruitcount.lf.mod = glmmTMB(fruit_count ~ local_foreign + mar_size_scaled 
                            + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                            family = nbinom2, 
                            ziformula = ~ 1,
                            data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
summary(fruitcount.lf.mod)
nrow(filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(fruit_count)))
# 3176

# seed number
seeds.lf.mod = glmmTMB(total_est_seeds ~ local_foreign + mar_size_scaled 
                       + (1|site/blk) + (1|sirepop/sire) + (1|dampop/dam), 
                       family = nbinom2, 
                       ziformula = ~ 1,
                      data = filter(plants_wi, mar_surv == 1, !is.na(mar_size_scaled), !is.na(total_est_seeds)))
summary(seeds.lf.mod)
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

results_la = results_new %>% 
  unite(model, model, part, remove = TRUE, sep = ".") %>% 
  gather(var, value, c(3:6)) %>% 
  unite(variable, var, rowname, sep = ":") %>% 
  spread(variable, value)

results_round = rapply(object = results_la, f = round, classes = "numeric", how = "replace", digits = 3) 

write.csv(results_round, "results/la_table.csv", row.names = FALSE)



# plotting predictions from overall model ---------------------------------

# raw mean
mean(plants_wi_seeds$total_est_seeds)
# 11.12115

# generate model predictions and CIs
pr.seedsall = ggpredict(seedsall.lf.mod, c("local_foreign"), ci.lvl = 0.95, typical = "mean")
pr.seedsall$x = c("foreign", "local")
pr.seedsall$predicted
# 11.51867 10.52239, reasonable match to raw mean

# generate raw averages per sire for plotting
real.seedsall = plants_wi %>% group_by(local_foreign, sire) %>% dplyr::summarize(mu = mean(total_est_seeds, na.rm = TRUE))

# plot predictions and raw sire averages
ggplot(data = pr.seedsall, aes(x = x, y = predicted)) +
  geom_jitter(data = real.seedsall, aes(x = local_foreign, y = mu), alpha = 0.4) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), size = 1.5) +
  ylab("Lifetime fitness") +
  xlab("Foreign vs. local origin")


