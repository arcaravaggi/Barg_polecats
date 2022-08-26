setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(lme4)
library(car)
library(ggplot2)
library(MuMIn)
library(tidyverse)
library(broom)
library(caret)
library(multcomp)
library(kableExtra)
library(data.table)
library(emmeans)

p.dat <- read.csv("polecat_model_data_lcm.csv")[,5:22]
str(p.dat)
head(p.dat)

# Read water data and stack
w.r <- read.csv("water_roadkill.csv")[,26]
w.p <- read.csv("water_pseudo.csv")[,17]

w.all <-  c(w.r, w.p)
w.all[is.na(w.all)] <- 0

# Convert to density of water features
p.dat$water_dens <- w.all/pi*1.5^2

# Rescale numeric variables between 0-1 so that they are directly comparable
scaleR <- function(x) (x-min(x))/(max(x) - min(x))
t <- data.frame(apply(p.dat[c(1:14,18:19)], 2, scaleR))

p.dat <- cbind(t, p.dat[,c(15:17)])

# drop observations where fclass = trunk_link
p.dat <- p.dat[!p.dat$fclass == "trunk_link",]

# Principal Components
p.pca <- preProcess(p.dat[,c(2:7, 9:10)], 
                       method=c("BoxCox", "center", "pca"),# Correct for skew 
                       thresh = 1) # retain components which cumulatively account for 100% of variance
p.comp <- predict(p.pca, p.dat[,c(2:7, 9:10)]) # components
p.load <- data.frame(p.pca$rotation) # loadings

# Add PCs to dataframe 
p.dat <- cbind(p.dat, p.comp)

# Extract PC variance
wdbc.pr <- prcomp(p.dat[,c(2:7, 9:10)], center = TRUE, scale = TRUE)
summary(wdbc.pr)

# Plot PCs
# Eigenvalue <1 explains less than single explanatory variable
screeplot(wdbc.pr, type = "l", npcs = 8, main = "Screeplot of the first 8 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

# Plot cumulative variance threshold
# Include 4 or 5?
cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
plot(cumpro[0:8], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(h = 0.75, col="blue", lty=5)
abline(v = 4, col = "blue", lty = 5)
legend("topleft", legend=c("Cut-off @ 90%"),
       col=c("blue"), lty=5, cex=0.6)

# Construct GLM
mod1 <- glm(presence ~ + elevation + patches + r_dens + water_dens + d_rabbit + PC1 + PC2 + PC3, 
            data = p.dat, 
            family = binomial(link = "logit"), na.action = na.fail)
summary(mod1)
vif(mod1) # Check Variance Inflation Factors (correlation)

# d_water is causing a separation problem. What rows are triggering the warning?
# glm0.resids <- augment(mod1) %>%
#  mutate(p = 1 / (1 + exp(-.fitted)),
#         warning=p>1-10 * .Machine$double.eps) %>%
#  subset(warning == T)
# They're where presence == 1. Clearly we can't drop these rows. 
# Drop d_water

# Construct average model
mod1.aic <- dredge(mod1, rank = "AIC")
mod1.topset <- subset(mod1.aic, delta <2) 
mod1.topmod <- subset(mod1.aic, delta == 0)

mod2 <- glm(presence ~ patches + elevation + r_dens + d_rabbit + PC1 + PC2, 
            data = p.dat, 
            family = "binomial"(link = logit), na.action = na.fail)
summary(mod2)
vif(mod2) # Check Variance Inflation Factors (correlation)


# Average model
mod1.avg <- model.avg(mod1.topset)
summary(mod1.avg)
importance(mod1.avg)

# Post-hoc between road classes
p.dat$fclass <- as.factor(p.dat$fclass)
mod3 <- glht(mod1, mcp(fclass = "Tukey"))
summary(glht(mod1, mcp(fclass = "Tukey")))
cld(mod3)
plot(mod3)

# Tables and figures ####

# Bar plot of deaths by roadtype
r.dat <- subset(p.dat, presence == 1)
r.dat <- data.frame(table(r.dat$fclass))
names(r.dat) <- c("class", "count")

r.dat$class <- str_to_title(r.dat$class)
r.dat$class[8] <- "Trunk link"

source('theme_ac1.R')
p1 <- ggplot(r.dat, aes(class, count)) + 
  geom_bar(stat = "identity") +
  xlab("Road Class") +
  ylab("Number of polecat mortalities") +
  theme_ac1(base_size_a = 14, base_size_t = 14) +
  theme(axis.text.x = element_text(angle = 25, vjust = 0.8, hjust=0.8)) +
  theme(text=element_text(family="sans"))

p1

ggsave(plot = p1, dpi = 300, filename = "../road_class.png", width = 6, height = 4, units = "in", device = "png")


# PC graphs
# Calculate mean and standard deviation
p.dat$Presence <- p.dat$presence
p.dat$Presence <- as.factor(p.dat$Presence)
pca.m <- p.dat %>% group_by(Presence) %>%
  summarise(avg.x = mean(PC1),
            avg.y = mean(PC2),
            sd.x = sd(PC1),
            sd.y = sd(PC2))

# Plot
p2 <- ggplot(pca.m, aes(avg.x, avg.y, colour = Presence)) +
  geom_point() +
  geom_errorbar(aes(ymin = avg.y-sd.y, ymax = avg.y+sd.y), width = 0.2) +
  geom_errorbar(aes(xmin = avg.x-sd.x, xmax = avg.x+sd.x), width = 0.2) +
  xlab("PC1 (20.5%)") +
  ylab("PC2 (17.1%)") +
  theme_ac1(base_size_a = 16, base_size_t = 16) + 
  theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size = 16))
  
ggsave(plot = p2, dpi = 300, filename = "../PC1_PC2.png", width = 6, height = 6, units = "in", device = "png")

# Make a pretty table of my pc loadings
p.load2 <- p.load
p.load2 <- (setattr(p.load2, "row.names", c("Grassland", "Wetland", "Arable",
                                            "Urban", "Conifer", "Broad-leaf", "Mixed",
                                            "Industrial")))
p.load2 <- rename(p.load2, "PC1 (21%)" = "PC1", "PC2 (17%)" = "PC2", "PC3 (14%)" = "PC3", "PC4 (13%)" = "PC4")

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

p.load2 <- round_df(p.load2, 3)


loadings <- p.load2[,1:4] %>% 
  kable() %>%
  kable_classic(full_width = T, html_font = "Calibri") %>%
  row_spec(0, bold = 1, align = "right")
  
loadings



# Traffic model ####
glm2 <- glm(presence ~ + elevation + patches + r_dens + all_mt_ + water_dens + d_rabbit + PC1 + PC2 + PC3, 
            data = p.dat, 
            family = binomial(link = "logit"), na.action = na.fail)
vif(glm2)
summary(glm2)

# Check for overdispersion. Ratio should be ~ 1
resid.ssq <- sum(residuals(glm2,type="pearson")^2)  ## sum of squares of Pearson resids
resid.df <- nrow(p.dat)-length(coef(glm2))        ## estimated resid df (N-p)
resid.ssq/resid.df    

# Get top subset and top model
mod2.aic <- dredge(glm2, rank = "AIC")
mod2.topset <- subset(mod2.aic, delta <2) 
mod2.topmod <- subset(mod2.aic, delta == 0)

top.mod <- glm(presence ~ + elevation + patches + r_dens + all_mt_ + d_rabbit + PC1 + PC2, 
               data = p.dat, 
               family = binomial(link = "logit"), na.action = na.fail)
summary(top.mod)
vif(mod2) # Check Variance Inflation Factors (correlation)

mod2.avg <- model.avg(mod2.topset)
summary(mod2.avg)
importance(mod2.avg)


