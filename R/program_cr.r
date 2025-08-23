## R function for censored linear regression

logit <- function(x) log(x)-log(1-x)

creg <- function(
  dat,
  tmb_file="censored_regression",
  half_life=c(2.1,30),       # 134Cs: 2.1 years, 137Cs: 30 years
  start_date=NULL,
  threshold=100,
  eval.max=1000,
  iter.max=1000,
  mis=0,
  ignore_censor=FALSE,
  est=TRUE,
  init_p=NULL,
  fix_sigma_re=FALSE,
  init_bb=c(5,5),
  init_cc=c(-6,6),
  init_dd=c(-8,-8),
  init_sigma_re=0.5,
  const_aa=c(1,2),     # constraint for log_Q
  const_bb=c(1,1),     # constraint for logit_P
  const_cc=c(1,1),     # constraint for log_k_1
  const_dd=c(1,1),     # constraint for log_k_2
  const_sigma=c(1,1),     # constraint for log_sigma
  one_component=FALSE,     # TRUE: one_component model, FALSE: two_component model
  silent=TRUE,
  do_compile=FALSE
){
require(TMB)

if (do_compile){
  compile(paste0(tmb_file,".cpp"))
  dyn.load(dynlib(tmb_file))
}

dat$Time <- as.Date(dat$Time)
dat <- dat[order(dat$Time),]

Cs_134 <- dat$Cesium134
Cs_137 <- dat$Cesium137
Censor_134 <- as.numeric(dat$Censor134)
Censor_137 <- as.numeric(dat$Censor137)
if (ignore_censor){
  Censor_134 <- Censor_137 <- rep(0, length(Censor_134))
}
Time <- ts(dat$Time)

start <- ifelse(is.null(start_date),Time[1],ts(as.Date(start_date)))

Time <- Time - start

n <- length(Cs_134) 

dat_for_est <- list(Cs_134=Cs_134, Cs_137=Cs_137, Censor_134=Censor_134, Censor_137=Censor_137, Time=as.numeric(Time), mis=mis, half_life=half_life, one_component=as.numeric(one_component))

if (one_component) {
  init_bb <- c(30,30)
  init_dd <- c(-30,-30)
  const_bb <- c(NA,NA)
  const_dd <- c(NA,NA)
}
HL <- log(0.5)/(half_life*365)

if (is.null(init_p)) {
  mod134 <- lm(log(Cs_134)~Time+offset(HL[1]*Time))
  mod137 <- lm(log(Cs_137)~Time+offset(HL[2]*Time))

  init_p[1:2] <- c(mod134$coef[1],mod137$coef[1])
  init_p[3:4] <- init_bb
  init_p[5:6] <- log(-c(mod134$coef[2],mod137$coef[2]))
  init_p[5:6] <- ifelse(is.na(init_p[5:6]),log(0.000001),init_p[5:6])
  init_p[7:8] <- init_dd
  init_p[9:10] <- c(log(summary.lm(mod134)$sigma),log(summary.lm(mod137)$sigma))
  init_p[11] <- log(init_sigma_re)
}

par_for_est <- list(
  aa = init_p[1:2],
  bb = init_p[3:4],
  cc = init_p[5:6],
  dd = init_p[7:8],
  log_sigma = init_p[9:10],
  log_sigma_RE = init_p[11],
  z = rep(0,n)
)

map_tmb <- NULL
map_tmb$aa <- factor(const_aa)
map_tmb$bb <- factor(const_bb)
map_tmb$cc <- factor(const_cc)
map_tmb$dd <- factor(const_dd)
map_tmb$log_sigma <- factor(const_sigma)
if (fix_sigma_re) {map_tmb$log_sigma_RE <- factor(NA); z <- rep(factor(NA),n)}

obj_tmb <- MakeADFun(dat_for_est, par_for_est, random=c("z"), map=map_tmb, DLL=tmb_file, silent=silent)

if (est) opt <- nlminb(obj_tmb$par, obj_tmb$fn, obj_tmb$gr, control=list(eval.max=eval.max, iter.max=iter.max)) else opt <- list(par=obj_tmb$par, objective=obj_tmb$fn(obj_tmb$par))

aa <- obj_tmb$report()$aa
bb <- obj_tmb$report()$bb
cc <- obj_tmb$report()$cc
dd <- obj_tmb$report()$dd
sigma_134 <- obj_tmb$report()$sigma[1]
sigma_137 <- obj_tmb$report()$sigma[2]
sigma_RE <- obj_tmb$report()$sigma_RE
Q <- obj_tmb$report()$Q
P <- obj_tmb$report()$P
K1 <- obj_tmb$report()$K1
K2 <- obj_tmb$report()$K2

parms <- c(aa, bb, cc, dd, sigma_134, sigma_137, sigma_RE)
names(parms) <- c("a_134", "a_137", "b_134", "b_137", "c_134", "c_137", "d_134", "d_137", "sigma_134", "sigma_137", "sigma_RE")

TwoCompPar <- rbind(Q, P, K1, K2)
rownames(TwoCompPar) <- c("Q", "P","K1","K2")
colnames(TwoCompPar) <- c("Cs134", "Cs137")

like <- -opt$objective
np <- length(opt$par)
aic <- c(like, np, (-2)*like+2*np)
names(aic) <- c("log_likelihood","n_pars","AIC")

list(obj=obj_tmb, opt=opt, start=start, mis=mis, ignore_censor=ignore_censor, half_life=half_life, one_component=one_component, fix_sigma_re=fix_sigma_re, const_aa=const_aa, const_bb=const_bb, const_cc=const_cc, const_dd=const_dd, const_sigma=const_sigma, parms=parms, TwoCompPar=TwoCompPar, aic=aic)
}

predict_cs <- function(
  res,
  x="2022/7/1",
  perc=0.95
){
  p <- res$parms
  Parms <- res$TwoCompPar
  
  Q <- Parms[1,]
  P <- Parms[2,]
  K1 <- Parms[3,]
  K2 <- Parms[4,]
  sigma_134 <- p["sigma_134"]
  sigma_137 <- p["sigma_137"]
  sigma_RE <- p["sigma_RE"]
  
  HL <- log(0.5)/(res$half_life*365)
 
  tx <- as.numeric(ts(as.Date(x)) - res$start)
 
  mean_log_cs_134 <- log(Q[1])+HL[1]*tx+log(P[1]*exp(-K1[1]*tx)+(1-P[1])*exp(-K2[1]*tx))
  mean_log_cs_137 <- log(Q[2])+HL[2]*tx+log(P[2]*exp(-K1[2]*tx)+(1-P[2])*exp(-K2[2]*tx))
  se_log_cs_134 <- sqrt(sigma_134^2+sigma_RE^2)
  se_log_cs_137 <- sqrt(sigma_137^2+sigma_RE^2)
  
  output <- rbind(
    c(mean_log_cs_134, se_log_cs_134, mean_log_cs_134-qnorm(1-(1-perc)/2)*se_log_cs_134, mean_log_cs_134+qnorm(1-(1-perc)/2)*se_log_cs_134),
    c(mean_log_cs_137, se_log_cs_137, mean_log_cs_137-qnorm(1-(1-perc)/2)*se_log_cs_137, mean_log_cs_137+qnorm(1-(1-perc)/2)*se_log_cs_137)
  )
  
  rownames(output) <- c("Cesium134", "Cesium137")
  colnames(output) <- c("mean","se","lo","up")
  
  output
}

# evaluation of risk

eval_risk <- function(
  res,
  tmb_file="risk_evaluation",
  threshold=100,
  Target_Date="2019/7/1",
  calc_se=TRUE,
  perc=0.95,
  eval.max=1000,
  iter.max=1000,
  simulate=FALSE,
  Sim=1000,
  B=1000,
  seed0=1,
  stats="mean",
  n_g=c(5000,150),
  lower_limit=10^(-6),
  fix_zero_risk=FALSE,
  silent=TRUE,
  do_compile=FALSE,
  dd=0.000001
){
require(stringr)
require(TMB)

if (do_compile){
  compile(paste0(tmb_file,".cpp"))
  dyn.load(dynlib(tmb_file))
}

start <- res$start
half_life <- res$half_life
one_component <- res$one_component

sdrep <- sdreport(res$obj)
parms <- sdrep$par.fixed
cov1 <- sdrep$cov.fixed

const_aa <- res$const_aa
const_bb <- res$const_bb
const_cc <- res$const_cc
const_dd <- res$const_dd
const_sigma <- res$const_sigma

aa_names <- c("a_134","a_137")
bb_names <- c("b_134","b_137")
cc_names <- c("c_134","c_137")
dd_names <- c("d_134","d_137")
sigma_names <- c("sigma_134","sigma_137")
log_sigma_names <- c("log_sigma_134","log_sigma_137")
par_names <- c(aa_names,bb_names,cc_names,dd_names,log_sigma_names,"log_sigma_RE") 

np <- length(parms)
nz <- length(par_names)

Alpha <- 1-perc
lmode <- NULL
    
set.seed(seed0)
if (fix_zero_risk) risk_est <- risk_se <- risk_lo <- risk_up <- 0 else{
  if (simulate){
    require(mvnfast)   
    
    if (calc_se){
      parms_z <- rmvn(B,parms,cov1)
      
      colnames(parms_z) <- names(parms)

      parms_zz <- matrix(NA, nrow=nrow(parms_z), ncol=nz)
      colnames(parms_zz) <- par_names
      parms_zz[,aa_names[which(is.na(const_aa))]] <- rep((res$parms[aa_names])[which(is.na(const_aa))],each=nrow(parms_zz))
      parms_zz[,aa_names[which(!is.na(const_aa))]] <- parms_z[,which(colnames(parms_z)=="aa")[na.omit(const_aa)]]
      parms_zz[,bb_names[which(is.na(const_bb))]] <- rep((res$parms[bb_names])[which(is.na(const_bb))],each=nrow(parms_zz))
      parms_zz[,bb_names[which(!is.na(const_bb))]] <- parms_z[,which(colnames(parms_z)=="bb")[na.omit(const_bb)]]
      parms_zz[,cc_names[which(is.na(const_cc))]] <- rep((res$parms[cc_names])[which(is.na(const_cc))],each=nrow(parms_zz))
      parms_zz[,cc_names[which(!is.na(const_cc))]] <- parms_z[,which(colnames(parms_z)=="cc")[na.omit(const_cc)]]
      parms_zz[,dd_names[which(is.na(const_dd))]] <- rep((res$parms[dd_names])[which(is.na(const_dd))],each=nrow(parms_zz))
      parms_zz[,dd_names[which(!is.na(const_dd))]] <- parms_z[,which(colnames(parms_z)=="dd")[na.omit(const_dd)]]
      parms_zz[,log_sigma_names[which(is.na(const_sigma))]] <- rep(log(res$parms[sigma_names])[which(is.na(const_sigma))],each=nrow(parms_zz))
      parms_zz[,log_sigma_names[which(!is.na(const_sigma))]] <- parms_z[,which(colnames(parms_z)=="log_sigma")[na.omit(const_sigma)]]
      if (!any(names(parms)=="log_sigma_RE")) parms_zz[,nz] <- log(res$obj$report()$sigma_RE) else parms_zz[,nz] <- parms_z[,ncol(parms_z)]
    
      z_s <- matrix(rnorm(B*Sim, 0, 1), nrow=B, ncol=Sim)
      w_134 <- matrix(rnorm(B*Sim, 0, 1), nrow=B, ncol=Sim)
      w_137 <- matrix(rnorm(B*Sim, 0, 1), nrow=B, ncol=Sim)
      sim_sample <- sapply(1:B, function(i) get(stats)(simulate_cs(res, Target_Date=Target_Date, Sim=Sim, output_option=3, par_set=parms_zz[i,], z_s=z_s[i,], w_134=w_134[i,], w_137=w_137[i,])[,2] > threshold))
  
      risk_est <- mean(sim_sample)
      risk_se <- sd(sim_sample)

      ci <- quantile(sim_sample, probs=c(Alpha/2,1-Alpha/2))
      risk_lo <- as.numeric(ci[1])
      risk_up <- as.numeric(ci[2])
    } else {
      parms_z <- parms
          
      parms_zz <- rep(NA, len=nz)
      names(parms_zz) <- par_names
      parms_zz[aa_names[which(is.na(const_aa))]] <- (res$parms[aa_names])[which(is.na(const_aa))]
      parms_zz[aa_names[which(!is.na(const_aa))]] <- parms_z[which(names(parms_z)=="aa")[na.omit(const_aa)]]
      parms_zz[bb_names[which(is.na(const_bb))]] <- (res$parms[bb_names])[which(is.na(const_bb))]
      parms_zz[bb_names[which(!is.na(const_bb))]] <- parms_z[which(names(parms_z)=="bb")[na.omit(const_bb)]]
      parms_zz[cc_names[which(is.na(const_cc))]] <- (res$parms[cc_names])[which(is.na(const_cc))]
      parms_zz[cc_names[which(!is.na(const_cc))]] <- parms_z[which(names(parms_z)=="cc")[na.omit(const_cc)]]
      parms_zz[dd_names[which(is.na(const_dd))]] <- (res$parms[dd_names])[which(is.na(const_dd))]
      parms_zz[dd_names[which(!is.na(const_dd))]] <- parms_z[which(names(parms_z)=="dd")[na.omit(const_dd)]]
      parms_zz[log_sigma_names[which(is.na(const_sigma))]] <- log(res$parms[sigma_names])[which(is.na(const_sigma))]
      parms_zz[log_sigma_names[which(!is.na(const_sigma))]] <- parms_z[which(names(parms_z)=="log_sigma")[na.omit(const_sigma)]]
      if (!any(names(parms)=="log_sigma_RE")) parms_zz[nz] <- log(res$obj$report()$sigma_RE) else parms_zz[nz] <- parms_z[length(parms_z)]
 
      z_s <- rnorm(Sim, 0, 1)
      w_134 <- rnorm(Sim, 0, 1)
      w_137 <- rnorm(Sim, 0, 1)
      risk_est <- mean(simulate_cs(res, Target_Date=Target_Date, Sim=Sim, output_option=3, par_set=parms_zz, z_s=z_s, w_134=w_134, w_137=w_137)[,2] > threshold)
      risk_se <- risk_lo <- risk_up <- NA
    }
  } else{
    require(statmod)
    if (length(n_g)==1) n_g <- rep(n_g, 2)
    gauss_hermite <- gauss.quad(n_g[1],kind="hermite")      # node and weight for gauss-hermite quadratue
    g_her <- cbind(gauss_hermite$nodes, gauss_hermite$weights)
    gauss_laguerre <- gauss.quad(n_g[2], "laguerre")      # node and weight for gauss-laguerre quadratue
    g_lag <- cbind(gauss_laguerre$nodes, exp(gauss_laguerre$nodes)*gauss_laguerre$weights)

    loc_aa <- which(names(parms) %in% "aa")
    loc_bb <- which(names(parms) %in% "bb")
    loc_cc <- which(names(parms) %in% "cc")
    loc_dd <- which(names(parms) %in% "dd")
    loc_sigma <- which(names(parms) %in% "log_sigma")

    target <- as.numeric(ts(as.Date(Target_Date)) - start)
    
    risk_f <- function(p,D,lmode=NULL,func_type=0,dd=rep(0,nz)){
      dat_for_risk <- list(node=g_lag[,1], wt=g_lag[,2], Time=target, D=D, func_type=func_type, half_life=half_life, one_component=as.numeric(one_component))

      par_for_risk <- list(
        aa=p[1:2]+dd[1:2],
        bb=p[3:4]+dd[3:4],
        cc=p[5:6]+dd[5:6],
        dd=p[7:8]+dd[7:8],
        log_sigma=log(p[9:10])+dd[9:10],
        log_sigma_RE=log(p[11])+dd[11],
        z=0
      )
      
      if (is.null(lmode)){
        obj_tmb <- MakeADFun(dat_for_risk, par_for_risk, map=map_tmb, DLL=tmb_file, silent=silent)
        opt <- nlminb(obj_tmb$par, obj_tmb$fn, obj_tmb$gr, control=list(eval.max=eval.max, iter.max=iter.max))
        sdrep_tmb <- sdreport(obj_tmb)
      
        muHat <- as.numeric(opt$par)
        sigmaHat <- as.numeric(sqrt(sdrep_tmb$cov.fixed))
      } else {
        muHat <- lmode$muHat
        sigmaHat <- lmode$sigmaHat
      }
      
      par_for_risk$z <- muHat

      obj_tmb <- MakeADFun(dat_for_risk, par_for_risk, map=map_tmb, DLL=tmb_file, silent=silent)
      
      prob_d <- 1-sqrt(2)*sigmaHat*sum(exp(log(g_her[,2])+g_her[,1]^2-sapply(muHat+sqrt(2)*sigmaHat*g_her[,1], obj_tmb$fn)))

      return(list(muHat=muHat, sigmaHat=sigmaHat, prob_d=prob_d))
    }
    
    map_tmb <- NULL
    map_tmb$aa <- rep(factor(NA),2)
    map_tmb$bb <- rep(factor(NA),2)
    map_tmb$cc <- rep(factor(NA),2)
    map_tmb$dd <- rep(factor(NA),2)
    map_tmb$log_sigma <- rep(factor(NA),2)
    map_tmb$log_sigma_RE <- factor(NA)
    
    risk_out <- risk_f(res$parms, threshold)
    risk_est <- risk_out$prob_d
    if (risk_est < lower_limit) risk_est <- 0
    lmode <- list(muHat=risk_out$muHat, sigmaHat=risk_out$sigmaHat)
  
    # derivatives
    
    if (calc_se){
      DD <- matrix(0, nrow=np, ncol=nz)
      
      loc_sigma_RE <- which(names(parms)=="log_sigma_RE")
      if (length(loc_sigma_RE) > 0) DD[loc_sigma_RE, nz] <- dd
      
      if (length(loc_aa)>0){
        r_pos <- which(names(parms)=="aa")
        c_pos <- which(str_detect(par_names,"a_"))
        nn <- length(r_pos)
        for (k in 1:nn) DD[r_pos[k],c_pos[which(const_aa==k)]] <- dd
      }
      if (length(loc_bb)>0){
        r_pos <- which(names(parms)=="bb")
        c_pos <- which(str_detect(par_names,"b_"))
        nn <- length(r_pos)
        for (k in 1:nn) DD[r_pos[k],c_pos[which(const_bb==k)]] <- dd
      }
      if (length(loc_cc)>0){
        r_pos <- which(names(parms)=="cc")
        c_pos <- which(str_detect(par_names,"c_"))
        nn <- length(r_pos)
        for (k in 1:nn) DD[r_pos[k],c_pos[which(const_cc==k)]] <- dd
      }
      if (length(loc_dd)>0){
        r_pos <- which(names(parms)=="dd")
        c_pos <- which(str_detect(par_names,"d_"))
        nn <- length(r_pos)
        for (k in 1:nn) DD[r_pos[k],c_pos[which(const_dd==k)]] <- dd
      }
      if (length(loc_sigma)>0){
        r_pos <- which(names(parms)=="log_sigma")
        c_pos <- which(str_detect(par_names,"sigma_1"))
        nn <- length(r_pos)
        for (k in 1:nn) DD[r_pos[k],c_pos[which(const_sigma==k)]] <- dd
      }
      
      risk_est_p <- risk_est_m <- NULL
      for (i in 1:np){
        risk_est_p <- c(risk_est_p, risk_f(res$parms, threshold, lmode=lmode, dd=+DD[i,])$prob_d)
        risk_est_m <- c(risk_est_m, risk_f(res$parms, threshold, lmode=lmode, dd=-DD[i,])$prob_d)
      }
      
      d_risk <- (risk_est_p-risk_est_m)/(2*dd)
    
      risk_se <- sqrt(t(d_risk)%*%cov1%*%d_risk)
    } else risk_se <- NA
      
    if (risk_est==0) {
        risk_lo <- risk_up <- 0
    } else{
        C1 <- exp(qnorm(1-Alpha/2)*risk_se/(risk_est*(1-risk_est)))

        risk_lo <- risk_est/(risk_est + (1-risk_est)*C1)
        risk_up <- risk_est/(risk_est + (1-risk_est)/C1)
    }
  }
}

risk <- data.frame(Prob=paste0("Prob(Cesium134+Cesium137 > ",threshold,") on ",Target_Date),Value=risk_est, Sim=simulate, SE=risk_se, Alpha=Alpha, Lower=risk_lo, Upper=risk_up)

return(risk)
}

plot_risk <- function(
res,
threshold=100,
start_date="2011/11/1",
end_date="2022/12/1",
len=10,
calc_se=TRUE,
perc=0.95,
simulate=FALSE,
Sim=1000,
B=1000,
seed0=1,
n_g=c(100,100),
lower_limit=10^(-6),
dd=0.00001,
continuity=FALSE,
log_scale=FALSE
){
require(tidyverse)

start_date <- as.Date(start_date)
end_date <- as.Date(end_date)

range_date <- seq(start_date, end_date, len=len)

risk_set <- sapply(range_date, function(x) eval_risk(res,threshold=threshold,Target_Date=x,calc_se=calc_se,perc=perc,simulate=simulate,Sim=Sim,B=B,seed0=seed0,lower_limit=lower_limit,dd=dd))

if (continuity){
  end <- which(sapply(1:len, function(i) unlist(risk_set[2,i])) < lower_limit)[1]
  risk_set[c(2,4,6:7),end:len] <- 0
}

y_label <- paste0("Prob(Cesium134+Cesium137 > ",threshold,")")

dat_for_plot <- data.frame(Date=range_date, Prob_low=as.numeric(risk_set[6,]), Prob_mid=as.numeric(risk_set[2,]), Prob_upp=as.numeric(risk_set[7,]))

p_tot <- ggplot(dat_for_plot,aes(x=Date,y=Prob_mid))+geom_ribbon(aes(ymin=Prob_low,ymax=Prob_upp),alpha=0.2)+geom_line(color="red")+labs(x="Date",y=y_label)+theme_bw()
if (log_scale) p_tot <- p_tot + scale_y_log10()

p_tot
}

plot_cs <- function(
dat,
res,
log_scale=TRUE,
perc=0.95,
len=100
){
require(tidyverse)

dat$Time <- as.Date(dat$Time)

tx <- as.Date(unique(dat$Time))

start_date <- min(tx)
end_date <- max(tx)

range_date <- seq(start_date, end_date, len=len)

pred <- sapply(range_date, function(x) predict_cs(res, x=x, perc=perc))

dat_cs134 <- data.frame(Time=range_date, mean=exp(pred[1,]), lo=exp(pred[5,]), up=exp(pred[7,]))
dat_cs137 <- data.frame(Time=range_date, mean=exp(pred[2,]), lo=exp(pred[6,]), up=exp(pred[8,]))

dat$Censor134 <- factor(dat$Censor134)
dat$Censor137 <- factor(dat$Censor137)

p_134 <- ggplot(dat_cs134,aes(x=Time,y=mean))+geom_ribbon(aes(ymin=lo,ymax=up),alpha=0.3)+geom_point(data=dat, aes(x=Time,y=Cesium134,color=Censor134),size=1)+geom_line(color="blue",size=1.5)+labs(x="Date",y="Cesium134")+theme_bw()
if (log_scale) p_134 <- p_134 + scale_y_log10()

p_137 <- ggplot(dat_cs137,aes(x=Time,y=mean))+geom_ribbon(aes(ymin=lo,ymax=up),alpha=0.3)+geom_point(data=dat, aes(x=Time,y=Cesium137,color=Censor137),size=1)+geom_line(color="blue",size=1.5)+labs(x="Date",y="Cesium137")+theme_bw()
if (log_scale) p_137 <- p_137 + scale_y_log10()

p_all <- cowplot::plot_grid(p_134, p_137)

p_all
}

simulate_cs <- function(
  res,
  Target_Date="2022/7/1",
  Sim=1000,
  output_option=1,
  par_set=NULL,
  z_s=NULL,
  w_134=NULL,
  w_137=NULL,
  seed0=1
){
# output_option: 1 = Cs_134 only, 2 = Cs_137 only, 3 = Total, 4 = All
require(TMB)

set.seed(seed0)

if (is.null(par_set)) parms_b <- log(res$parms) else parms_b <- par_set

HL <- log(0.5)/(res$half_life*365)

a_134 <- parms_b[1]
a_137 <- parms_b[2]
b_134 <- parms_b[3]
b_137 <- parms_b[4]
c_134 <- parms_b[5]
c_137 <- parms_b[6]
d_134 <- parms_b[7]
d_137 <- parms_b[8]
sigma_134 <- exp(parms_b[9])
sigma_137 <- exp(parms_b[10])
sigma_RE <- exp(parms_b[11])
p_134 <- 1/(1+exp(-b_134))
p_137 <- 1/(1+exp(-b_137))

start <- as.numeric(res$start)

target <- as.numeric(ts(as.Date(Target_Date)) - start)

if (is.null(z_s)) z_s <- rnorm(Sim, 0, 1)
z <- sigma_RE*z_s

if (output_option!=2) {
  if (is.null(w_134)) w_134 <- rnorm(Sim, 0, 1)
  cs_134 <- z+a_134+HL[1]*target+log(p_134*exp(-exp(c_134)*target)+(1-p_134)*exp(-exp(d_134)*target))+sigma_134*w_134
} else cs_134 <- NA
if (output_option!=1) {
  if (is.null(w_137)) w_137 <- rnorm(Sim, 0, 1)
  cs_137 <- z+a_137+HL[2]*target+log(p_137*exp(-exp(c_137)*target)+(1-p_137)*exp(-exp(d_137)*target))+sigma_137*w_137
} else cs_137 <- NA

Cs_134 <- exp(cs_134)
Cs_137 <- exp(cs_137)

if (output_option==1) output <- data.frame(z=z, Cs_134=Cs_134)
if (output_option==2) output <- data.frame(z=z, Cs_137=Cs_137)
if (output_option==3) output <- data.frame(z=z, Cs_total=Cs_134+Cs_137)
if (output_option==4) output <- data.frame(z=z, Cs_134=Cs_134, Cs_137=Cs_137, Cs_total=Cs_134+Cs_137)

return(output)
}

risk_date_search <- function(
  p,
  res,
  threshold=100, 
  n_g=c(100,100), 
  calc_se=FALSE,
  init_log_Time=2,
  min_Date="2011/3/11"
){
  risk <- function(x) (eval_risk(res, threshold=threshold, Target_Date=as.Date(min_Date)+exp(x), n_g=n_g, calc_se=FALSE)$Value-p)^2

  opt <- nlminb(init_log_Time, risk)

  Time <- as.Date(min_Date)+exp(opt$par)
  
  risk <- eval_risk(res, threshold=threshold, Target_Date=Time, n_g=n_g, calc_se=calc_se)
  
  list(opt=opt,min_Date=as.Date(min_Date),est_Date=Time,risk=risk)
}

plot_cs2 <- function(
  dat,
  res,
  res1,
  res2,
  log_scale=TRUE,
  perc=0.95,
  len=100
){
  dat_fp <- dat
  dat_fp$Time <- as.Date(dat_fp$Time)

  tx <- unique(dat_fp$Time)

  start_date <- as.Date(tx[1])
  end_date <- as.Date(tx[length(tx)])

  range_date <- seq(start_date, end_date, len=len)

  pred <- sapply(range_date, function(x) predict_cs(res, x=x, perc=perc))
  pred1 <- sapply(range_date, function(x) predict_cs(res1, x=x, perc=perc))
  pred2 <- sapply(range_date, function(x) predict_cs(res2, x=x, perc=perc))

  dat_cs134 <- data.frame(Time=range_date, mean=exp(pred[1,]), lo=exp(pred[5,]), up=exp(pred[7,]))
  dat_cs137 <- data.frame(Time=range_date, mean=exp(pred[2,]), lo=exp(pred[6,]), up=exp(pred[8,]))

  dat_fp$Censor134 <- ifelse(dat_fp$Censor134==0, "No", "Yes")
  dat_fp$Censor137 <- ifelse(dat_fp$Censor137==0, "No", "Yes")

  dat_fp$Censor134 <- factor(dat_fp$Censor134)
  dat_fp$Censor137 <- factor(dat_fp$Censor137)

  p_134 <- ggplot(dat_cs134,aes(x=Time,y=mean))+geom_ribbon(aes(ymin=lo,ymax=up),alpha=0.1)+geom_point(data=dat_fp, aes(x=Time,y=Cesium134,color=Censor134, shape=Censor134),size=1.5)+geom_line(color="blue",size=1.5)+geom_line(aes(x=Time,y=exp(pred1[1,])),color="purple3",size=1.5,linetype="dashed")+geom_line(aes(x=Time,y=exp(pred2[1,])),color="brown",size=1.5,linetype="dotted")+labs(x="Date",y="Cesium134", color="Censored")+theme_bw()
if (log_scale) p_134 <- p_134 + scale_y_log10() + theme(legend.position = "none")

  p_137 <- ggplot(dat_cs137,aes(x=Time,y=mean))+geom_ribbon(aes(ymin=lo,ymax=up),alpha=0.1)+geom_point(data=dat_fp, aes(x=Time,y=Cesium137,color=Censor137, shape=Censor137),size=1.5)+geom_line(color="blue",size=1.5)+geom_line(aes(x=Time,y=exp(pred1[2,])),color="purple3",size=1.5,linetype="dashed")+geom_line(aes(x=Time,y=exp(pred2[2,])),color="brown",size=1.5,linetype="dotted")+labs(x="Date",y="Cesium137", color="Censored")+theme_bw()
if (log_scale) p_137 <- p_137 + scale_y_log10() + theme(legend.position = "none")

  cowplot::plot_grid(p_134, p_137, ncol=2)
}

halflife_estimator <- function(
  res,
  p0=c(200,200),
  dd=10^(-6)
){
  one_component <- res$one_component
  sdrep <- sdreport(res$obj)
  if (one_component){
    p <- sdrep$par.fixed[3]
    cov_fixed <- sdrep$cov.fixed[3,3]  
  } else{
    p <- sdrep$par.fixed[c(4:5,3)]
    cov_fixed <- sdrep$cov.fixed[c(4:5,3),c(4:5,3)]
  }
   
  HL <- log(0.5)/(res$half_life*365)
 
  obj <- function(log_tx, par, HL){
    tx <- exp(log_tx)
    K <- exp(par)
    
    mean_log_cs <- HL*tx-K*tx
     
    (mean_log_cs-log(0.5))^2
  }
  
  obj_a <- function(log_tx, par, HL){
    tx <- exp(log_tx)
    K1 <- exp(par[1])
    K2 <- exp(par[2])
    P <- 1/(1+exp(-par[3]))
    
    mean_log_cs <- HL[1]*tx+log(P*exp(-K1*tx)+(1-P)*exp(-K2*tx))
     
    (mean_log_cs-log(0.5))^2
  }
  
  DD <- diag(length(p))
  diag(DD) <- dd
  
  Out <- NULL
  for (k in 1:2){
    Est <- NULL
    for (i in 1:min(2,length(p))){
      output <- nlm(obj, log(p0[k]), par=p[i], HL=HL[k])
      output_P <- nlm(obj, output$estimate, par=p[i]+dd, HL=HL[k])
      output_N <- nlm(obj, output$estimate, par=p[i]-dd, HL=HL[k])
      D_out <- (exp(output_P$estimate)-exp(output_N$estimate))/(2*dd)
      if (one_component) SE <- abs(D_out)*sqrt(cov_fixed) else SE <- abs(D_out)*sqrt(cov_fixed[i,i])
      Est <- rbind(Est, c(exp(output$estimate), SE))
    }
    
    if (!one_component){
      output <- nlm(obj_a, log(p0[k]), par=p, HL=HL[k])
      output_P <- sapply(1:length(p), function(j) nlm(obj_a, output$estimate, par=p+DD[j,], HL=HL[k])$estimate)
      output_N <- sapply(1:length(p), function(j) nlm(obj_a, output$estimate, par=p-DD[j,], HL=HL[k])$estimate)
      D_out <- (exp(output_P)-exp(output_N))/(2*dd)
      SE <- sqrt(D_out%*%cov_fixed%*%D_out)
      Est <- rbind(Est, c(exp(output$estimate), SE))
    }
    
    Out <- rbind(Out, Est)
  }
  
  if (one_component) rownames(Out) <- c("Cs134","Cs137") else  rownames(Out) <- c("Cs134:Fast","Cs134:Slow","Cs134:Total","Cs137:Fast","Cs137:Slow","Cs137:Total")
  colnames(Out) <- c("Estimate","SE")
  
  Out
}
