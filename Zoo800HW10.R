#########Zoology800 Homework 10 - Maddy Teter - 11/6/2025

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)

#Objectiggplot2#Objective 1 - Use own data and select two variables that share a relationship with more than 
#30 paired datapoints

#Read in csv files (comparing agb to sediment grain size)
aboveground_biomass <- read.csv("PSD25_aboveground_biomass.csv")
sediment_grain_size <- read.csv("PSD25_sediment_grain_size.csv")
#Load sediment BD data (need this to get grain size proportions)
sediment_bulk_density <- read.csv("PSD25_sediment_bulk_density.csv")

#Se function to apply to values
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length((x)))

#calculate final agb
aboveground_biomass <- aboveground_biomass %>%
  mutate(final_dry_mass_g = dry_mass_g - dish_mass_g, na.rm = TRUE)

#filter for thalassia only at a single time point
aboveground_biomass_average <- aboveground_biomass%>%
  filter(species == c("Tt"), time_interval_wks == "t7")

#Calculate sed BD final mass, mean and se
sediment_bulk_density <- sediment_bulk_density %>%
  mutate(final_bd_gml = (dry_mass_g - wp_mass_g)/sample_vol_ml)

#Calculate mean and se for sedbd (mg/L)
sediment_bulk_density_average <- sediment_bulk_density %>%
  group_by(site_id) %>%
  summarise(
    mean_sedbd_gml = mean(final_bd_gml, na.rm = TRUE),
    se_sedbd_gml = se(final_bd_gml))

#Create column for final dw (g) for sediment grain sizes
sediment_grain_size <- sediment_grain_size %>%
  mutate(final_sedgs_dw_g = dry_mass_g - wp_mass_g)

#Join mean_sedbd column from sedbd df to sedgs df
sediment_grain_size <- sediment_grain_size %>%
  left_join(sediment_bulk_density_average %>% select(site_id, mean_sedbd_gml, 
                                 se_sedbd_gml), by = "site_id")


#Pivot to wide format — one row per sample_id/site_id with gravel type dw listed
#as multiple columns
sediment_grain_size_wide <- sediment_grain_size %>%
  pivot_wider(
    id_cols = c(sample_id, mean_sedbd_gml),  # keep sample-specific columns
    names_from = grain_type,
    values_from = final_sedgs_dw_g,
    names_prefix = "",
    values_fn = list(final_sedgs_dw_g = mean)  # in case of duplicates
  )

#Calculate fine sediment grain size
#((mean_sedbd_gml*60(sedgs sample size)) - gravel_dw_g - sand_dw_g = fine_dw_g
#for each sample id/ grain size sample collected
sediment_grain_size_wide <- sediment_grain_size_wide %>%
  mutate(
    fine = ifelse(
      is.na(fine),
      (mean_sedbd_gml * 60) - gravel - sand,
      fine
    )
  )

#Join back to original data
sediment_grain_size_updated <- sediment_grain_size %>%
  pivot_wider(names_from = grain_type, values_from = final_sedgs_dw_g) %>%
  left_join(sediment_grain_size_wide, by = c("sample_id", "mean_sedbd_gml")) %>%
  mutate(
    gravel = coalesce(gravel.y, gravel.x),
    sand = coalesce(sand.y, sand.x),
    fine = coalesce(fine.y, fine.x)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

#Change back to long format (grain size type is one column)
sediment_grain_size_final <- sediment_grain_size_updated %>%
  pivot_longer(
    cols = c(gravel, sand, fine),
    names_to = "grain_type",
    values_to = "final_sedgs_dw_g"
  )

#Create column for proportion_sedgs_dw for each grain type in each sample
sediment_grain_size_final <- sediment_grain_size_final %>%
  mutate(proportion_sedgs_dw_g = final_sedgs_dw_g/(mean_sedbd_gml*60))

#Filter for fine sedgs only
fine_sediment_grain_size <- sediment_grain_size_final %>%
  filter(grain_type == c("fine"))

#Join fine sed gs proportion to aboveground biomass df
aboveground_biomass_average <- aboveground_biomass_average %>%
  left_join(fine_sediment_grain_size %>% select(site_id, proportion_sedgs_dw_g),
            by = "site_id")

#Collapse df to one row per site_id*column replicate
aboveground_biomass_average <- aboveground_biomass_average %>%
  group_by(site_id, replicate) %>%  # or whatever your replicate column is named
  summarize(
    across(
      everything(),
      ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else first(.x)),
    .groups = "drop")

#Run linear model to compare relationship between agb and fine sed grain size proportion
aboveground_biomass_fine <- lm(final_dry_mass_g ~ proportion_sedgs_dw_g, 
                               data = aboveground_biomass_average)
summary(aboveground_biomass_fine)
#Results - significant relationship (p-value: 0.00148)

#Plot relationship - it is linear
ggplot(aboveground_biomass_average, aes(x = final_dry_mass_g, y = proportion_sedgs_dw_g)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  labs(
    title = "Relationship between Thalassia Aboveground Biomass and Fine Grain Sediment Proportion",
    x = "Aboveground Biomass (g)",
    y = "Fine Sediment Grain Size Proportion"
  ) +
  theme_classic()

#Check residuals
autoplot(aboveground_biomass_fine)
#Homoscedastic residuals - they don't funnel or change vertical distance from 0 line by much
#along x-axis. There is a dip around 0.9 and a slight rise around 1.5 but nothing drastic.
#This indicates that the variance is mostly equal 
#There are a few outliers (11, 33, 35) that are highlighted in the plot and could influence
#the data slightly

#check for normality
hist(aboveground_biomass_fine$residuals)
#normal distribution

#No autocorrelation - each point is independent from each other 

#Create median and 95% quartile predictions
median <- median(aboveground_biomass_average$proportion_sedgs_dw_g, na.rm = TRUE)
quartile_95 <- quantile(aboveground_biomass_average$proportion_sedgs_dw_g, 0.95, na.rm = TRUE)

#Create new df 
aboveground_biomass_fine_predictions <- data.frame(proportion_sedgs_dw_g = 
                                                     c(median, quartile_95))

#predict(object, newdata, interval = “prediction”)
predict(aboveground_biomass_fine, aboveground_biomass_fine_predictions, 
        interval = "prediction")

#The model has a lot of variance in general, but the results show that the 95% quartile range
#has a larger difference between lower and upper values compared to the median which reflects
#less certainty in the values. Median predicts lower mass average (0.95) compared to 
#95% quartile prediction (1.32).


#Objective 2
set.seed(99) #uses same set of values each time the model is ran 

#Set paramters for linear model
n <- 100
beta_0 <- 1.0 #intercept
beta_1 <- 2.0 #slope

#Generate predictor x values
x <- runif(n, 0, 10)

#Create lognormal error values (unimodal but not normal)
error <- rlnorm(n, meanlog = 0, sdlog = 0.4) - exp(0 + 0.4^2 / 2)

#establish linear model
y <- beta_0 + beta_1 * x + error

#fit the linear model
linear_model_1 <- lm(y ~ x)
summary(linear_model_1)


#Repeat the process to compare true and estimated slope and intercept
set.seed(99)
#New parameter 
n_simulations <- 100 #run simulations 100 times

simulation_results <- replicate(n_simulations, { #use same code as before but nested with n_simulations
  x <- runif(n, 0, 10)
  error <- rlnorm(n, meanlog = 0, sdlog = 0.4) - exp(0 + 0.4^2 / 2)
  y <- beta_0 + beta_1 * x + error
  linear_model_2 <- lm(y ~ x)
  coefs <- coef(linear_model_2)
  
#Generate predictor x values
simulation_predictions <- predict(linear_model_2, newdata = data.frame(x = x),
                   interval = "prediction", level = 0.95)
  
#find the % coverage of how often the prediction values fall in 95% CI
  coverage <- mean(y >= simulation_predictions[, "lwr"] & 
                     y <= simulation_predictions[, "upr"])
  
#return results from entire simulation (coef 1 = estimated intercept, coef 2 = 
#estimated slope)
  c(
    est_intercept = coefs[1],
    est_slope = coefs[2],
    coverage = coverage
  )
}, simplify = TRUE)

#Create dataframe with all information from estimated values
final_simulation_results <- as.data.frame(t(simulation_results))

#The intercept and slope values from the estimated model stay fairly close to the true values
#of 1 (intercept) and 2 (slope). 84 out of the 100 estimated values fall inside of 95% 
#coverage - so majority of values are withing CI.This means that the estimated uncertainty
#is pretty close to the true uncertainty estimate. 



