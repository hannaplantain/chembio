install.packages("ggpubr")
install.packages("MASS")    
install.packages("kSamples")
remove.packages("Matrix")
utils::install.packages("lme4", type = "source")
utils::install.packages("Matrix", type = "source")
library(ggpubr)
library(MASS)    
library(lme4)
library(Matrix)
library(dplyr)
library(ggplot2)

CONC_I2M2 <- read.csv("CONC_I2M2.csv")
TU_I2M2_sum <- read.csv("TU_I2M2_sum.csv")
long_CONC_TU_I2M2 <- read.csv("long_CONC_TU_I2M2.csv")
long_Sum_TU <- read.csv("sum_TU.csv")
rm(Sum_TU)
DF1 <- TU_I2M2_sum

#add id 
DF1$id <- seq_along(DF1$Sample_Source_Code)
DF1 <- DF1[, c("id", setdiff(names(DF1), "id"))]

DF2 <- long_Sum_TU
DF2$id <- seq_along(DF2$Sample_Source_Code)
DF2 <- DF2[, c("id", setdiff(names(DF2), "id"))]
DF2 <- DF2 %>%
  mutate(id = as.character(id))

#scatter
DF2_clean <- DF2 %>%
  filter(is.finite(Sum_TU), is.finite(I2M2))

ggplot(DF_outlierfree, aes(x = Sum_TU, y = I2M2)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method = "loess", formula = y ~ x) + 
  theme_minimal() +
  labs(
    title = "Toxic pressure and I2M2",
    x = "Sum_TU",
    y = "I2M2"
  )
install.packages("rmcorr")
library(rmcorr)

rmc_ti <- rmcorr(Sample_Source_Code, Sum_TU, I2M2, DF_outlierfree)
print(rmc_ti)

DF_outlierfree <- DF_outlierfree %>%
  filter(!is.na(Sum_TU), !is.na(I2M2))

rmc_plot <- ggplot(DF_outlierfree,
                   aes(x = Sum_TU, y = I2M2, color = factor(Sample_Source_Code))) +
  geom_rmc(rmc_ti) +
  geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, color = "black")  +
  theme_minimal() +
  theme(legend.position="none") +
  labs(
    title = "Repeated measure correlation between toxic pressure and I2M2",
    x = "Sum_TU",
    y = "I2M2"
  )

print(rmc_plot)
plot(rmc_ti, overall = TRUE, lty = 2)
plot(rmc_ti)

?
#removing outliers

##case by case
DF2 <- DF2 %>%
  filter(id != "6902")
DF2 <- DF2 %>%
  filter(id != "8604")
DF2 <- DF2 %>%
  filter(id != "768")

str(DF2)
remove_outliers <- function(data) {
  numeric_columns <- data %>% select(where(is.numeric))  # Only numeric columns
  
  # Apply outlier removal to numeric columns
  numeric_columns_outlierfree <- numeric_columns %>%
    mutate(across(everything(), ~ {
      Q1 <- quantile(., 0.25, na.rm = TRUE)
      Q3 <- quantile(., 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      # Replace outliers with NA to keep vector length
      ifelse(. >= lower_bound & . <= upper_bound, ., NA)
    }))
  
  # Combine with non-numeric columns
  non_numeric_columns <- data %>% select(where(~ !is.numeric(.)))
  data_no_outliers <- bind_cols(non_numeric_columns, numeric_columns_outlierfree)
  
  return(data_no_outliers)
}
DF_outlierfree <- remove_outliers(DF2)

rlanghead(DF1)
str(DF1)

#again, scatterplot
ggplot(DF1, aes(x = `330-54-1`, y = I2M2)) + 
  geom_point() + 
  theme_minimal()

#explore normality of data - result - I2M2 not normally distributed
ggdensity(DF1$I2M2,  
          main = "density I2M2",
          xlab = "I2M2") 

ggqqplot(DF$I2M2) #light-tailed both sides
shapiro.test(DF1$I2M2) #significant - not normal


ggplot(DF1, aes(x = Sum_TU, y = I2M2)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE, color = "purple") + 
  theme_minimal()

#correlation -----
cor(DF$Sum_TU, DF$I2M2, method = "spearman") # one year - preselect

##aggregate
DF_aggregated <- DF1 %>%
  group_by(Sample_Source_Code, Year) %>%
  summarise(
    Sum_TU = mean(Sum_TU, na.rm = TRUE),  
    I2M2 = mean(I2M2, na.rm = TRUE),
    .groups = "drop"
  )
##correlation between Sum_TU and I2M2
correlation <- cor(DF_aggregated$Sum_TU, DF_aggregated$I2M2, use = "complete.obs")
print(correlation)

#model
lm <- lmer(I2M2  ~ Sum_TU + Year + (1 | Sample_Source_Code), data = long_Sum_TU)
summary(lm)
lm
