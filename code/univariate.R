######################################################################
## R code generating graphs and numbers of the first two sections of
## the book chapter "Prospective Detection of Outbreaks" by
## B. Allévius and M. Höhle in the Handbook of Infectious Disease Data Analysis.
##
## Author: Michael Höhle <http://www.math.su.se/~hoehle>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 2017-11-10
######################################################################

##Load packages
library(ggplot2)
library(tidyverse)
library(magrittr)
library(surveillance)
library(lubridate)

######################################################################
## Create plot of time series of monthly number of meningococcal cases
## in Sect. 24.1
######################################################################

##Load data pre-processed from the surveillance package
meningo <- readRDS("meningo_monthly.rds")

##Aggregate over space
meningo_ts <- meningo %>% group_by(year, month) %>%
  summarise(count=sum(count), population=sum(population),
            time=min(time),date=min(date)) %>%
  mutate(quarter=factor(paste0("Q",quarter(date)))) %>%
  ungroup


##Make a plot over time
ts <- ggplot( meningo_ts, aes(x=date, y=count)) +
  geom_crossbar(aes(ymin=0,ymax=count),fatten=0,fill="lightgray") +
  xlab("Time (month)") + ylab("No. cases") + ylim(c(0,NA)) +
  scale_x_date(date_breaks = "6 month",
               date_minor_breaks = "1 month",
               labels = scales::date_format("%b-%Y")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1))

##Look at it
ts

##Store as pdf on disk
ggsave(file.path("..","ChapterInBook","chapters","chapter23","figures","meningo-monthly-ts.pdf"), ts, width = 8, height = 4)

######################################################################
## Text about significance levels of Ears and Farrington methods in
## Sect. 24.2
######################################################################

##Number of observations
k <- 7
##True alpha when considering t-distribution and prediction factor s.t.
##the multiplication factor for sd is 3.
alpha <- 1- pt(3/sqrt(1+1/k), df=k-1)
alpha
sprintf("%.4f",alpha)
##Sanity check. Result should be 3.
qt(1-alpha, df=k-1) * sqrt(1+1/k)

##alpha when using multiplication factor of 3 and assuming underlying Gaussian:
alpha <- 1-pnorm(3)
sprintf("%.4f",alpha)
##If we really want this alpha, then the multiplication factor should be:
sprintf("%.2f",qt(1-alpha, df=k-1) * sqrt(1+1/k))

######################################################################
## Farrington plot in Fig. 24.2 which corresponds to the prediction
## interval for Jan 2008
######################################################################

##Add info on what to monitor and what is historic values
meningo_ts %<>% mutate(monitor = lubridate::year(date) >  2007)

the_t <- min(meningo_ts %>% filter(monitor) %$% time)
the_date <- min(meningo_ts %>% filter(monitor) %$% date)

##Define the set of time points to include as historic values
historic_times <- the_t - rep((1:3)*12, each=2*3+1) + rep((-3):3,3)

##Prepare glm data (todo: take this directly form meningo)
glm_data <- meningo_ts %>% filter(time %in% c(the_t,historic_times)) %>%
  mutate(historic= time %in% historic_times, pi_lower=0,pi_upper=0)

##Fit glm to thie historic data
m <- glm(count ~ 1 + time, data=glm_data %>% filter(historic),
         family=quasipoisson)

##Trend is only included, if significant?
summary.glm(m)$coefficients["time","Pr(>|t|)"]

##Extract overdispersion parameter
phi.hat <- summary(m)$dispersion

##Make a data frame with new data for visualization
newdata <- data.frame(time = seq(min(glm_data$time),the_t,by=0.01)) %>%
  mutate(date=seq(min(glm_data$date),the_date,length=n()),
         historic = date <= max(glm_data %>% filter(historic) %$% date))

##Note the choice of alpha! (we pick the small value which one assumes
##by considering the mean, i.e. by not knowing the k dependent
##prediction factor as well as the use of the t-distribution alpha
alpha <- 1-pnorm(3)
alpha

##Make the predictions
p <- predict(m, newdata=newdata, type="response",se=TRUE)
newdata %<>% mutate(fit = p$fit,
                    se = sqrt( phi.hat * fit + p$se.fit^2),
                    pi_lower = fit - qnorm(1-alpha) * se,
                    pi_upper = fit + qnorm(1-alpha) * se)

##Information about the monitored time point
newdata %>% filter(date == the_date)

##Upper limit of the prediction interval
ut <- sprintf("%.1f",newdata %>% filter(date == the_date) %$% pi_upper)
##Actual observation
yt <- glm_data %>% filter(date == the_date) %$% count
##Sound an alarm?
data.frame(ut=ut, yt=yt, alarm=yt>ut)

##Illustrate everything in a plot
p71 <- ggplot(newdata %>% filter(date == the_date), aes(x=date,ymin=pi_lower,ymax=pi_upper,color=historic)) + geom_errorbar(width=60) +
  theme_bw() +
  geom_point(data=glm_data, aes(x=date, y=count,color=historic,shape=historic)) +
  geom_line(data=newdata,aes(x=date,y=fit,color=historic,linetype=historic)) +

  xlab("Time (month)") + ylab("No. cases") +
  ##scale_color_grey(name = "Historic?") +
  scale_color_manual(name = "Historic?", values=c("gray5","gray65")) +
  scale_linetype_manual(name = "Historic?", values=c(3,1,6)) +
  scale_shape_discrete(name = "Historic?")


##Show plot
p71

##Store to disk
ggsave(file.path("..","ChapterInBook","chapters","chapter23","figures","farrington-predict71.pdf"), p71, width = 6, height = 2.75)

######################################################################
## Same analysis but now using the surveillance package
######################################################################

sts <- with(meningo_ts, sts(observed=count, epoch=as.numeric(date), frequency=12,epochAsDate=TRUE))
plot(sts)

##Compute and extract the upper limit of the prediction interval
s_farrington <- farrington(sts, control=list(range=the_t,alpha=(1-pnorm(3))*2,b=3,w=3, reweight=FALSE,powertrans="none"))
##This should be 21.1 but is 21.7, because a trend is not included!
sprintf("%.1f",upperbound(s_farrington))
##Value of the upper limit when using a 2/3-power transformation
s_farrington_23 <- farrington(sts, control=list(range=the_t,alpha=(1-pnorm(3))*2,b=3,w=3, reweight=FALSE,powertrans="2/3"))
sprintf("%.1f",upperbound(s_farrington_23))

##value of the negbin quantile
mu <- newdata %>% filter(date == the_date) %$% fit
qnbinom(1-alpha, mu=mu, size=1/((phi.hat-1)/mu))
qnbinom(1-alpha, mu/(phi.hat-1), 1/phi.hat)



######################################################################
##Compare for b=1 with EARS
######################################################################

y <- meningo_ts %>% slice(the_t -12 + ((-3):3)) %$% count
##7 obs when w=3 and b=3
#y <- unlist(meningo_ts[min(which(meningo_ts$monitor))-rep(12*(1:3),each=7) + rep((-3):3,times=3),"count"])
y
sprintf("%.1f",mean(y))
sprintf("%.1f",sd(y))

##Select alpha
alpha <- 1-pnorm(3)
##alpha, if we respect that it's a prediction interval
##alpha <- 1- pt(3/sqrt(1+1/k), df=k-1)

##Limits
sprintf("%.1f",mean(y) + 3*sd(y))
k <- length(y)
sprintf("%.1f",mean(y) + qt(1-alpha, df=k-1) * sqrt(1+1/k) * sd(y))

##Farrington outputs - no transformation
s_farrington <- farrington(sts, control=list(range=the_t,alpha=alpha*2,b=1,w=3, reweight=FALSE,powertrans="none"))
sprintf("%.1f",upperbound(s_farrington))
##value of the 2/3 power transformation limit
s_farrington_23 <- farrington(sts, control=list(range=the_t,alpha=alpha*2,b=1,w=3, reweight=FALSE,powertrans="2/3"))
sprintf("%.1f",upperbound(s_farrington_23))
s_farrington_23 <- farrington(sts, control=list(range=the_t,alpha=alpha*2,b=1,w=3, trend=FALSE,reweight=FALSE,powertrans="2/3"))
sprintf("%.1f",upperbound(s_farrington_23))

##Fit the glm to see overdispersion parameter is <1.
b <- 1
historic_times <- the_t - rep((1:b)*12, each=2*3+1) + rep((-3):3,b)
glm_data <- meningo_ts %>% filter(time %in% c(the_t,historic_times)) %>%
  mutate(historic= time %in% historic_times, pi_lower=0,pi_upper=0)

##No time trend, because there are not 3 years of data
m <- glm(count ~ 1, data=glm_data %>% filter(historic), family=quasipoisson)
(phi.hat <- max(1,summary(m)$dispersion))
##Fit as poisson since dispersion param is 1
m <- glm(count ~ 1, data=glm_data %>% filter(historic), family=poisson)

##Fit with negative binomial doesn't work.
if (FALSE) {
  require(mgcv)
  m_nb <- gam(count ~1, family=nb(theta=NULL), data=glm_data %>% filter(historic))
  summary(m_nb)
  m_nb <- MASS::glm.nb(count ~ 1, data=glm_data %>% filter(historic))
  summary(m_nb)
}

##Prediction
p <- predict(m, newdata=glm_data %>% filter(!historic),type="response",se=TRUE)

##The prediction
(mu <- p$fit)
p$fit + qnorm(1-alpha) * sqrt( phi.hat * p$fit + p$se.fit^2)

threshold <- if (phi.hat<=1) {
               qpois(1-alpha,mu)
             } else {
               qnbinom(1-alpha, mu/(phi.hat-1), 1/phi.hat)
             }
threshold

qpois(1-alpha,mu)


##Specificity, i.e. P(no alarm | X from an in-control distribution)
##  P( X <= threshold | X from an in-control distribution)
##So the FPF is the false alarm rate, i.e.
##  P( X > threshold | X from an in-control distribution) =
##  1 - P( X <= threshold | X from an in-control distribution)
sprintf("%.3f%%",alpha*100)
sprintf("%.3f",(1-ppois(mean(y) + 3*sd(y),mu))*100)
1-ppois(mean(y) + qt(1-alpha, df=k-1) * sqrt(1+1/k) * sd(y), mu)
1-ppois(upperbound(s_farrington_23),mu)
1-ppois(qpois(1-alpha,mu),mu)
alpha

##What if observations are only a quarter of their current size?
meningo_ts2 <- meningo_ts %>% mutate(count= round(count/4))
y2 <- meningo_ts2 %>% slice(the_t -12 + ((-3):3)) %$% count
##False alarm probability of the EARS approach
sprintf("%.2f%%",(1-ppois(mean(y2) + 3*sd(y2), mean(y2)))*100)
##Of the farrington with quantile threshold
1-ppois(qpois(1-alpha, mean(y2)), mean(y2))
