# packages

library(TMB)
library(tidyverse)
library(statmod)
library(mvnfast)

# data import

dat <- read.csv("example_dat.csv")
dat <- dat %>% mutate(across(Time, as.Date))

# install program

source("program_cr.r")

# activate TMB programs

tmb_file <- "censored_regression"
compile(paste0(tmb_file,".cpp"))
dyn.load(dynlib(tmb_file))

tmb_file <- "risk_evaluation"
compile(paste0(tmb_file,".cpp"))
dyn.load(dynlib(tmb_file))

# run program

res11 <- creg(dat, mis=0, one_component=TRUE, const_sigma=c(1,1))   # one-component model handling censored data via cumulative distribution functions
res12 <- creg(dat, mis=1, one_component=TRUE, const_sigma=c(1,1))   # one-component model ignoring censored data
res13 <- creg(dat, ignore_censor=TRUE, one_component=TRUE, const_sigma=c(1,1))    # one-component model treating censored data as uncensored.

res21 <- creg(dat, const_sigma=c(1,1))   # two-component model handling censored data via cumulative distribution functions
res22 <- creg(dat, mis=1, const_sigma=c(1,1))   # two-component model ignoring censored data
res23 <- creg(dat, ignore_censor=TRUE, const_sigma=c(1,1))    # two-component model treating censored data as uncensored.
  
# evaluate risk

D <- c(10,50,100)
ND <- length(D)
TD <- as.Date(c("2017/04/01"))
n_g <- c(1500,150)
risk21_a <- sapply(1:ND, function(i) eval_risk(res21, threshold=D[i], Target_Date=TD, calc_se=TRUE, n_g=n_g, simulate=FALSE))
risk21_s <- sapply(1:ND, function(i) eval_risk(res21, threshold=D[i], Target_Date=TD, calc_se=TRUE, Sim=1000, B=1000, simulate=TRUE))

# plot cesium contaminations

min_T <- min(dat$Time)
max_T <- max(dat$Time)
Times <- seq(min_T, max_T, len=500)
pred_a <- sapply(1:length(Times), function(i) predict_cs(res21,Times[i])[,1])

dat2 <- data.frame(Time=Times, pred_134=exp(pred_a[1,]), pred_137=exp(pred_a[2,]))

p11 <- ggplot(dat, aes(x=Time, y=Cesium134))+geom_point(size=2,aes(shape=Censor134, color=Censor134))+geom_line(data=dat2, aes(x=Time, y=pred_134), linewidth=1.5, color="blue")+scale_y_log10(labels = scales::label_number(accuracy=0.1))+labs(title=bquote(.("a) ") ~.("Example")~ .(": ")~{{}^{134}*Cs}),y=expression(paste({{}^{134}},"Cs",sep="")))+theme_classic()+theme(legend.position = "none")

p12 <- ggplot(dat, aes(x=Time, y=Cesium137))+geom_point(size=2,aes(shape=Censor137, color=Censor137))+geom_line(data=dat2, aes(x=Time, y=pred_137), linewidth=1.5, color="blue")+scale_y_log10(labels = scales::label_number(accuracy=0.1))+labs(title=bquote(.("b) ") ~.("Example")~ .(": ")~{{}^{137}*Cs}),y=expression(paste({{}^{137}},"Cs",sep="")), shape="Censored", color="Censored")+theme_classic()+theme(legend.position = c(0.9,0.9), legend.justification = c(1,1))

p1 <- cowplot::plot_grid(p11, p12)

# plot risk projection

p_risk1_1 <- plot_risk(res21, threshold=D[1], start_date="2011/3/11", end_date=TD, simulate=TRUE, len=50)+labs(y=bquote(.("Prob(")~{{}^{134}*Cs}+{{}^{137}*Cs}~.("> ")~.(D[1])~.("Bq/kg)")))
p_risk1_2 <- plot_risk(res21, threshold=D[2], start_date="2011/3/11", end_date=TD, simulate=TRUE, len=50)+labs(y=bquote(.("Prob(")~{{}^{134}*Cs}+{{}^{137}*Cs}~.("> ")~.(D[2])~.("Bq/kg)")))
p_risk1_3 <- plot_risk(res21, threshold=D[3], start_date="2011/3/11", end_date=TD, simulate=TRUE, len=50)+labs(y=bquote(.("Prob(")~{{}^{134}*Cs}+{{}^{137}*Cs}~.("> ")~.(D[3])~.("Bq/kg)")))

p2 <- cowplot::plot_grid(p_risk1_1, p_risk1_2, p_risk1_3, nrow=1)
