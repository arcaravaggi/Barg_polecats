# Library 
library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)
library(mgcViz)

# working directory
setwd("~/Dropbox/Polecat data/CSVdata")

# data import
freq <- read.csv("../CSVdata/seasonal_roadkill.csv")
head(freq)

# format dates
freq$date <- as.Date(freq$date, format="%d/%m/%Y")


# Do some models
# First one without a correlation argumnet
m1 <- gamm(obs ~ s(Month, bs = "cc", k = 12) + s(total_months), family=poisson,
          data = freq)
summary(m1$gam) #r2 = 0.471, F = 2.787, t = 22.87


# Do some more models but with a correlation argument for the month variable
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
m2 <- gamm(obs ~ s(Month, bs = "cc", k = 12) + s(total_months, k = 20),
           data = freq, correlation = corARMA(form = ~ 1|Year, p = 1), family = poisson,
           control = ctrl)
summary(m2$gam) #r2 = 0.469

m3 <- gamm(obs ~ s(Month, bs = "cc", k = 12) + s(total_months, k = 20),
           data = freq, correlation = corARMA(form = ~ 1|Year, p = 2), family = poisson,
           control = ctrl)
summary(m3$gam) #R2 = 0.466

m4 <- gamm(obs ~ s(Month, bs = "cc", k = 12) + s(total_months, k = 20),
           data = freq, correlation = corARMA(form = ~1|Year, p = 3), family = poisson,
           control = ctrl)
summary(m4$gam) # R2 = 0.461

m5 <- gamm(obs ~ s(Month, bs = "cc", k = 12) + s(total_months, k = 20),
           data = freq, correlation = corARMA(form = ~1|Year, p = 4), family = poisson,
           control = ctrl)
summary(m5$gam) # R2 = 0.433

# compare models
AIC(m1$lme) 
AIC(m2$lme)
AIC(m3$lme)
AIC(m4$lme)
AIC(m5$lme)

# Plot models
layout(matrix(1:2, ncol = 2))
plot(m1$gam, scale = 0)
layout(1)

layout(matrix(1:2, ncol = 2))
res <- resid(m1$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
pacf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1)

want <- seq(1, nrow(freq), length.out = 24)
pdat <- with(freq,
             data.frame(total_months = total_months[want], date = date[want],
                        Month = Month[want]))

## predict trend contributions
p1 <- predict(m1$gam, newdata = pdat, type = "terms", se.fit = TRUE)
p2 <- predict(m2$gam, newdata = pdat, type = "terms", se.fit = TRUE)
p3 <- predict(m3$gam, newdata = pdat, type = "terms", se.fit = TRUE)
p4 <- predict(m4$gam, newdata = pdat, type = "terms", se.fit = TRUE)

## combine with the predictions data, including fitted and SEs
pdat <- transform(pdat,
                  p1 = p1$fit[,2], se1 = p1$se.fit[,2],
                  p2 = p2$fit[,2], se2 = p2$se.fit[,2],
                  p3 = p3$fit[,2], se3 = p3$se.fit[,2],
                  p4 = p4$fit[,2], se4 = p4$se.fit[,2])

op <- par(mar = c(5,4,2,2) + 0.1)
ylim <- with(pdat, range(p1, p2, p3, p4))
ylim[1] <- floor(ylim[1])
ylim[2] <- ceiling(ylim[2])

plot(obs - mean(obs) ~ date, data = freq, type = "n",
     ylab = "Number of roadkill", ylim = ylim, xlab="Date")
lines(p1 ~ date, data = pdat, col = "red")
lines(p2 ~ date, data = pdat, col = "blue")
lines(p3 ~ date, data = pdat, col = "forestgreen", lwd = 1)
lines(p4 ~ date, data = pdat, col = "yellow") # Wow v curvy line, very bad.

legend("topleft",
       legend = c("Uncorrelated Errors", paste0("AR(", 1:4, ") Errors")),
       bty = "n", col = c("black","red","blue","forestgreen", "yellow"),
       lty = 1, lwd = c(1,1,1))
par(op)


layout(matrix(1:2, ncol = 2))
plot(m1$gam, scale = 0)
layout(1)

# make a pretty plot of the best model (m1)
b <- getViz(m1$gam)
o <- plot(sm(b, 1))
o + l_fitLine(colour = "black") + #l_rug(mapping = aes(x=x, y=y), alpha = 5.0) +
  l_ciLine(level = 0.95, mul = 5, colour = "black", linetype = 2) +
  scale_x_continuous(name="Month", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12), 
                     labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                     limits=c(1, 12))+ geom_hline(yintercept = 0, lty=3) +
  ylab("Smoothed number of detections") + 
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(colour = "black", size=1))


