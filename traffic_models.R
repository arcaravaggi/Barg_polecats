setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(multcomp)
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


# Data handling ####
# Mortality, pseudoabsence and traffic data were joined with roads
mort <- read.csv("./CSVdata/roads_roadkill.csv")
pseu <- read.csv("./CSVdata/roads_pseudo.csv")
traff <- read.csv("../CSVdata/roads_traffic.csv")

mort2 <- data.frame(presence = 1, # Extract relevant traffic and road data
                     mort[,11:17])

pseu2 <- data.frame(presence = 0,# Extract relevant traffic and road data
                    pseu[,3:9])

dat <- rbind(mort2, pseu2) # Combine mortality and pseudoabsence data

# Merge combined data with traffic join using longitude and latitude columns that were carried over from specific roads
dat2 <- merge(dat, traff, by = c("join_longi", "join_latit"))
dat2 <- dat2[,c(3,9,13,15)]
names(dat2) <- c("presence", "year", "traffic_vol", "road_class")
dat2$road_class <- as.factor(dat2$road_class)

# See how many records of each type remain
table(dat2$presence)

# Simple GLM exploring traffic volume and road class
glm1 <- glm(presence ~ traffic_vol + road_class, data = dat2)
summary(glht(glm1, mcp(road_class="Tukey")))


# Data from other modelling procedure #### 
# These data need to be combined with the traffic data to support the full modelling procedure
p.dat <- read.csv("../CSVdata/polecat_model_data.csv")[,5:21]
str(p.dat)

# Read water data and stack
w.r <- read.csv("../CSVdata/water_roadkill.csv")[,26]
w.p <- read.csv("../CSVdata/water_pseudo.csv")[,17]

w.all <-  c(w.r, w.p)
w.all[is.na(w.all)] <- 0

# Convert to density of water features
p.dat$water_dens <- w.all/pi*1.5^2

dat <- rbind(mort2, pseu2) # Combine mortality and pseudoabsence data

p.dat2 <- cbind(p.dat, dat) # Combine original model data with traffic data from above

p.dat2 <- p.dat2[,-c(14, 15, 19, 20, 22, 26)]

# Merge combined data with traffic join using longitude and latitude columns that were carried over from specific roads
p.dat3 <- merge(p.dat2, traff, by = c("join_longi", "join_latit"))
p.dat3 <- p.dat3[,-c(1, 2, 19:24, 26, 28:30)] # Remove unnecessary columns
colnames(p.dat3)[c(17,18)] <- c("traffic_vol", "road_class")

# See how many records of each type remain
table(p.dat3$presence)

# Rescale numeric variables between 0-1 so that they are directly comparable
scaleR <- function(x) (x-min(x))/(max(x) - min(x))
t <- data.frame(apply(p.dat3[c(1:8, 10:13, 16, 17)], 2, scaleR))
names(t) <- names(p.dat3[,c(1:8, 10:13, 16, 17) ] )

p.dat4 <- cbind(p.dat3[,c(9, 14, 15, 18)], t)

# Principal Components
p.pca <- preProcess(p.dat4[,c(5:12)], 
                    method=c("BoxCox", "center", "pca"),# Correct for skew 
                    thresh = 1) # retain components which cumulatively account for 100% of variance
p.comp <- predict(p.pca, p.dat4[,c(5:12)]) # components
p.load <- data.frame(p.pca$rotation) # loadings

# Add PCs to dataframe and drop original habitat covariates
p.dat4 <- cbind(p.dat4, p.comp)
hab.vars <- p.dat[,c(5:12)]
p.dat4[,c(5:12)] <- NULL


# GLM of presence by all variables, minus Moran's I
glm2 <- glm(presence ~ . - moran - d_road - d_building - d_water - PC5 - PC6 - PC7 - PC8, data = p.dat4, family = binomial(), na.action = 'na.fail')
vif(glm2)
summary(glm2)

# Check for overdispersion. Ratio should be ~ 1
resid.ssq <- sum(residuals(glm2,type="pearson")^2)  ## sum of squares of Pearson resids
resid.df <- nrow(p.dat4)-length(coef(glm2))        ## estimated resid df (N-p)
resid.ssq/resid.df    

# Get top subset and top model
mod2.aic <- dredge(glm2, rank = "AIC")
mod2.topset <- subset(mod2.aic, delta <2) 
mod2.topmod <- subset(mod2.aic, delta == 0)

top.mod <- glm(presence ~ patches + PC1 + PC2 + PC4 + road_class + elevation + traffic_vol,
            data = p.dat4, 
            family = "binomial"(link = logit), 
            control = list(maxit = 50), na.action = na.fail)
summary(top.mod)
vif(mod2) # Check Variance Inflation Factors (correlation)

mod2.avg <- model.avg(mod2.topset)
summary(mod2.avg)
importance(mod2.avg)

