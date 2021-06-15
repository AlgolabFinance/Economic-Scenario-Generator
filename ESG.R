rm(list=ls())
dev.off(dev.list()["RStudioGD"])

#   term : Term of the zero coupon bond of data input (1m,1y,2y,3y,5y,7y,10y,12y,15y,20y,30y,40y,50y,60y)
#   r    : The interest rate used to generate the next interest rate.
#   r0   : The initial interest rate. 
#   kappa: The mean reversion rate. 
#   theta: The mean rate or long term rate. 
#   sigma: Volatility. 
#   dt   : The change in time between observations. Defaults to 1/252 because
#          we assume generation of daily rates and there are 252 trading days 
#          per year. 
#   N    : The number of points to generate in each simulation. For example, 
#          the number of days when simulating daily rates.
#   M    : The number of simulations to run. 
#   years: The length or maturity of the bond.  
#   delta: Increment of maturity for the partial derivative calculation.
#   input_date: Date of calibration for Hull White Extended Vasicek model
#   data_range: Range of historical data uses for Vasicek parameters

library(dplyr)
library(ggplot2)
library(BBmisc) #for is.error
setwd("D:/OneDrive/1.Allianz/1.ESG/Code")

#__________________Define model parameters___________________
term <- 1 #EUR ZCB rate observed at t=0 maturity 1 year (1m,1y,2y,3y,4y,5y,7y,10y,12y,15y,20y,30y,40y,50y,60y)
input_date <-"31/12/2019"
years <- 10
N <- years * 252 # each year consists of 252 days
t <- (1:N)/252 # for plotting purposes
max.maturity=60
dt = 1/252
M <- 5 # test with several M simultions
set.seed(666)
delta <- 0.001
date_range = 133 #start date at 31/12/2014, data_origin start at 30/06/2014

#__________________________________Functions________________________________
VasicekShortRate <- function(mat) {
  # Simulate a vector of short rates using the Vasicek model.
  N <-mat*252
  short.rates <- rep(0, N)
  short.rates[1] <- r0
  for (i in 2:N) {
    # Calculates the next rate based on the discretization of the HullWhite model.
    short.rates[i] <- (short.rates[i - 1] + kappa*theta*dt + sigma*rnorm(n=1))/(1+kappa*dt)
  }
  return(short.rates)
}

HWShortRate1y <- function(theta,kappa,sigma) {
  # Simulate a vector of short rates at time t=1y using the Hullwhite extended Vasicek model
  N <-252
  short.rates <- rep(0, N)
  short.rates[1] <- r00
  for (i in 2:N) {
    # Calculates the next rate based on the discretization of the HullWhite model.
    short.rates[i] <- (short.rates[i - 1] + kappa*theta*dt + sigma*rnorm(n=1))/(1+kappa*dt)
  }
  return(short.rates[252])
}

Simulations <- function(N,M) {
  # Generates several short rate simulations 
  # Return an N row by M column matrix of simulated short rates. 
  sim.matrix <- matrix(nrow = N, ncol = M)
  for (i in 1:M) {
    sim.matrix[, i] <- VasicekShortRate(N/252)
  }
  return(sim.matrix)
}

VasicekZCBPrice <- function(mat) {
  # Calculates th zero coupon bond price at observed at time t=1 maturity T-t. 
  # using Vasicek model
  B.vas <- (1-exp(-mat*kappa)) / kappa 
  A.vas <- exp((B.vas-mat)*(kappa^2*theta-(sigma^2/2))/(kappa^2)-(sigma^2)/(4*kappa)*B.vas^2)
  ZCBPrize <- A.vas*exp(-B.vas*VasicekShortRate(mat)[mat*252]/100)
  return(ZCBPrize)
}

VasicekYieldCurve <- function() {
  # Produces a yield curve from the Hullwhite model with maturities ranging 
  # from 1 year to max.maturity.  
  # Returns:
  #   A decimal price of the bond (i.e. 0.98 for 98). 
  yields <- rep(0, 15)
  prices <- rep(0, 15)
  for (mat in 1:15) {
    prices[mat] <- VasicekZCBPrice(maturity[mat,1])
    yields[mat] <- -log(VasicekZCBPrice(maturity[mat,1]))/maturity[mat,1]
  }
  result <- data.frame(yields,prices)
  return(result)
}

VasicekParameters <- function(term) {
  # Calibrates the Vasicek model using the maximum likelihood estimator
  # Returns a vector of the form c(kappa, theta, sigma, r0), where kappa is the mean reversion rate, theta is the long-term rate/mean, sigma is the volatility and r0 is the last observed rate.
  # Later we shall use these parameter for HullWhite extended Vasicek model
  # which means the mean reversion rate and sigma are constants (time-independent) as Vasicek model
  data <- data.frame (data_all[term])
  n <- (nrow(data))-1
  S0 <- sum(data[1:n,1])/n
  S1 <- sum(data[2:(n+1),1])/n
  S00 <- as.numeric(crossprod(data[1:n,1], 
                              data[1:n,1]))/n
  S01 <- as.numeric(crossprod(data[1:n,1], 
                              data[2:(n+1),1]))/n
  S11 <- as.numeric(crossprod(data[2:(n+1),1], 
                              data[2:(n+1),1]))/n
  theta  <- (S1 * S00 - S0 * S01) / ( S0*S1 - S0^2 - S01 + S00)
  kappa <- log((S0 - theta) / (S1 - theta )) / dt
  beta <- (1-exp(-kappa*dt))/kappa
  sigma2 <- (S11-2*S01+S00-2*kappa*beta*(theta*S1-S01-theta*S0+S00)+
               kappa^2*beta^2*(n*theta^2-2*theta*S0+S00))/(n*beta*(1-kappa*beta/2))
  sigma <- sqrt(sigma2)
  r0 <- as.numeric(data[input_date,])
  
  return(c(r0, theta, kappa, sigma))
}

#__________________input data and PCA___________________
data_origin <- read.csv("Data.csv", sep=";",header=TRUE, row.name="Date")
data_origin <- na.omit(data_origin)
r00 <- as.numeric(data_origin[input_date,1])
data_origin <- data_origin[(date_range):nrow(data_origin),]

# Apply PCA to yield curve
PCA_matrix <- data.frame((data_origin %>% select(1:15)))
#data_all <- PCA_matrix
V <- princomp(PCA_matrix)$loadings[,c("Comp.1","Comp.2","Comp.3")]
data_all <- data.frame((as.matrix(PCA_matrix) %*% V) %*% t(V))
maturity <- data.frame("maturity" = c(0.083333,	1,	2,	3,	4,	5,	7,	10,	12,	15,	20,	30,	40,	50,	60))
summary (princomp(PCA_matrix))

ggplot(data.frame(V),aes(x=maturity[,1],y=Comp.1))+geom_line(color="red")+
  geom_line(data.frame(V),mapping=aes(x=maturity[,1],y=Comp.2),color="blue")+
  geom_line(data.frame(V),mapping=aes(x=maturity[,1],y=Comp.3),color="green")+
  geom_point(data.frame(V),mapping=aes(x=maturity[,1],y=Comp.1),color="red",size=3)+
  geom_point(data.frame(V),mapping=aes(x=maturity[,1],y=Comp.2),color="blue",size=2)+
  geom_point(data.frame(V),mapping=aes(x=maturity[,1],y=Comp.3),color="green",size=1.5)+
  ggtitle('First 3 principal components')+
  theme(plot.title = element_text(hjust = 0.5, vjust=2.12, lineheight=.8, face="bold")) +
  xlab("maturity") + # modif boot
  ylab("")  

#____________________
#Calculate the instantaneous forward rate see at time 0 maturity T
# as the partial derivative of ZCB price at time 0 maturity T with respect to T
calcul <- data.frame(t(data_origin[input_date,] %>% select(1:15)))
calcul <- data.frame(calcul,maturity)
names(calcul) <-c("rate","maturity")
calcul$priceZCB <- exp(-calcul$rate * calcul$maturity/100)

# Find all Vasicek parameter
para_table <-data.frame()
for (i in 1:15){
  para <- VasicekParameters(i) 
  para_table <- rbind(para_table,data.frame(t(para)))
}

names(para_table) <- c( "r0", "theta", "kappa", "sigma")
para_table$t <- calcul$maturity

#_________Plot 1 path and simulation of M paths using Vasicek parameter_______________

kappa <- para_table$kappa[term]
sigma <- para_table$sigma[term]
r0 <- para_table$r0[term]
theta <- para_table$theta[term]
data <- data.frame (data_all[term],rownames(data_all))
names(data) <- c("rate","date")
data$date <- as.Date(as.character(row.names(data)), "%d/%m/%Y")
start_date <- as.Date(as.character(input_date), "%d/%m/%Y")

matrix <-data.frame(Simulations(N,M))
projection_date <- as.numeric(rownames(matrix))+start_date
matrix$projection_date <-projection_date

test <- data.frame("short_rate"=VasicekShortRate(years),projection_date)
ggplot(test,aes(x=projection_date,y=short_rate))+  
  geom_line(alpha=0.8)+
  ggtitle('Short rate example')+     
  xlab("date") + 
  ylab("rate")

projection <- data.frame()
for (i in 1:(length(matrix)-1)){
  matrixi <- data.frame("scenario"= i,"rate"=matrix[,i],"date"=projection_date)
  projection <-rbind(projection,matrixi) 
}

quantiles <- data.frame(do.call("rbind", tapply(projection$rate, 
    projection$date, quantile, probs = c(0.10,0.50,0.90))))

p2<-ggplot(projection,aes(x=date,y=rate))+
  geom_line(data,mapping=aes(x=date,y=rate), color="red", size = 1) +
  ggtitle('HullWhiteprojection')+
  geom_line(aes(x=date,y=rate,group=scenario,color = as.factor(scenario)),size=0.3) +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5, vjust=2.12, lineheight=.3, face="bold")) +
  xlab("date de projection") + 
  ylab("taux 1 an - EUR") +ylim(-1,5)
p2+ geom_smooth(quantiles,mapping=aes(x=projection_date,y=quantiles[,3]), linetype = "dashed", color="brown", size = 1)+
  geom_smooth(quantiles,mapping=aes(x=projection_date,y=quantiles[,2]), linetype = "dashed", color="blue", size = 1)+
  geom_smooth(quantiles,mapping=aes(x=projection_date,y=quantiles[,1]), linetype = "dashed", color="green", size = 1)

#__________________Hull-White one factor model_________________
mon_spline <- function(x)
{
  return (t(spline(c(1/12.,1:5,7,10,12,15,20,30,40,50,60),x,xout=0:120)$y))
}
cubic_curve <- data.frame("rate"=t(mon_spline(calcul$rate)))
cubic_curve$maturity=c(0:120)
cubic_curve$maturity[1] <- 1/12

cubic_curve$price <- exp(-cubic_curve$rate*cubic_curve$maturity/100)
plot(cubic_curve$maturity,cubic_curve$price)
calcul3 <- data.frame(cubic_curve$maturity[1:61],
                      cubic_curve$maturity[2:62],
                      cubic_curve$price[1:61],
                      cubic_curve$price[2:62])
names(calcul3) <-c("t","t1","p_t","p_t1")
calcul3$price_minusdelta <- (calcul3$p_t1 - delta *
                               (calcul3$p_t1 -calcul3$p_t)/
                               (calcul3$t1-calcul3$t))
calcul3$price_plusdelta <- (calcul3$p_t+ delta *
                              (calcul3$p_t1 -calcul3$p_t)/
                              (calcul3$t1-calcul3$t))

Instfwdrate = -(log(calcul3$price_plusdelta[1])-log(cubic_curve$price[1]))/delta
calcul2 <- data.frame(-(log(calcul3$price_plusdelta[2:61])-log(calcul3$price_minusdelta[1:60]))/(2*delta))
calcul3$Instfwdrate <- data.frame(rbind(c(Instfwdrate),calcul2))
names(calcul3[,7]) <-c("Instfwdrate2")

cubic_price <-maturity
cubic_price$x <- data.frame("current"=calcul$priceZCB[1:15])
cubic_price[1,1] <-0
for (i in 3:62){
  cubic_price[i] <- data.frame(c(rep(0,15)))
}
rownames(cubic_price) <-c("0","1y","2y","3y","4y","5y","7y","10y","12y","15y","20y","30y","40y","50y","60y")
names(cubic_price) <- c("mat",c(0:60))
ESG_HW <- cubic_price[,3:62]
ESG_HW_1m <- ESG_HW[1,]
for (i in 3:62){
  for (j in 1:15){
    cubic_price[j,i] <- cubic_curve$price[i+cubic_price$mat[j]-1]
  }
}

Vasicek_ShortRate <- VasicekShortRate(60)
para_table$vasicekB <-(1-exp(-para_table$kappa*(para_table$t)))/para_table$kappa

for (i in 1:60){
  short_rate <- Vasicek_ShortRate[i*252]
  ESG_HW[,i] <- cubic_price[,i+2]/cubic_price[1,i+2]*
    exp(para_table$vasicekB*calcul3$Instfwdrate[i+1]-
          1/(4*para_table$kappa)*(para_table$sigma^2)*
          (1-exp(-2*para_table$kappa*(i)))*(para_table$vasicekB^2))*
    exp(-para_table$vasicekB*short_rate/100)
  ESG_HW[1,i] <- exp(-short_rate/12/100)
}

#__________________________Derive a yield curve_________________________
ESG_Vasicek <- data.frame(t(read.csv("EC_44_10000_RN_PRIIPs-simus_avg.csv",header=TRUE)[2,3:17]))
ESG_Vasicek$rate <- -log(ESG_Vasicek[,1])/maturity*100
names(ESG_Vasicek) <- c("PCA+BootstrapYC")

VasicekYC <- data.frame(VasicekYieldCurve())

Yieldcurves <- data.frame("maturity"=maturity,
                          "CurrentYC"=as.numeric(data_origin[input_date,]%>%select(1:15)),
                          "HullWhiteYC"=-log(ESG_HW$`1`)/maturity*100,
                          "VasicekYC" = VasicekYC[1],
                          "PCABootstrapYC"=ESG_Vasicek[,2])
names(Yieldcurves) <- c("maturity","CurrentYC","HullWhiteYC","VasicekYC","PCABootstrapYC")
plot(Yieldcurves$maturity, VasicekYC[,1], xlab = 'Maturity', type = 'l', ylab = 'Yield',
     main = 'VasicekYC',ylim = c(min(VasicekYC[,1]), max(VasicekYC[,1])),
     col = 1) 
plot(Yieldcurves$maturity, Yieldcurves$PCABootstrapYC,, xlab = 'Maturity', type = 'l', ylab = 'Yield',
     main = 'BootstrapYC',
     col = 1) 
plot(Yieldcurves$maturity, Yieldcurves$CurrentYC, xlab = 'Maturity', type = 'l', ylab = 'Yield',
     main = 'currentYC',
     ylim = c(min(Yieldcurves$CurrentYC), max(Yieldcurves$CurrentYC)),col = 1)
plot(Yieldcurves$maturity, Yieldcurves$HullWhiteYC, xlab = 'Maturity', type = 'l', ylab = 'Yield',
     main = '',
     ylim = c(min(Yieldcurves$HullWhiteYC), max(Yieldcurves$HullWhiteYC)),col = 1) 


#___________________________Validation - market consistency_________________
Yieldcurves$BondPrice_Current <- exp(-Yieldcurves$CurrentYC*Yieldcurves$maturity/100)
Yieldcurves$BondPrice_HW <- ESG_HW$`1`
Yieldcurves$BondPrice_Vasicek <- exp(-Yieldcurves$VasicekYC*Yieldcurves$maturity/100)
Yieldcurves$BondPrice_PCABootstrap <- exp(-Yieldcurves$PCABootstrapYC*Yieldcurves$maturity/100)
RMSE_HW <- sqrt(sum((Yieldcurves$BondPrice_HW- Yieldcurves$BondPrice_Current)^2)/nrow(Yieldcurves))
RMSE_Vasicek <- sqrt(sum((Yieldcurves$BondPrice_Vasicek- Yieldcurves$BondPrice_Current)^2)/nrow(Yieldcurves))
RMSE_PCABootstrap <- sqrt(sum((Yieldcurves$BondPrice_PCABootstrap- Yieldcurves$BondPrice_Current)^2)/nrow(Yieldcurves))

ggplot(Yieldcurves,aes(x=maturity,y=CurrentYC))+
  geom_line(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_HW), color="blue", size = 0.7) +
  geom_line(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_PCABootstrap), color="green", size = 0.7) +
  geom_line(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_Current), color="black", size = 0.7) +
  geom_point(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_HW), color="blue", size = 1.2) +
  geom_point(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_PCABootstrap), color="green", size = 1.2) +
  geom_point(Yieldcurves,mapping=aes(x=maturity,y=BondPrice_Current), color="black", size = 1.2) +
  ggtitle('')+
  xlab("Maturity") + 
  geom_line(deflators,mapping=aes(x=Yieldcurves$maturity,y=expected_price), color="red", size = 0.7) +
  ylab("ZCB Price")

print(c(RMSE_HW,RMSE_PCABootstrap, RMSE_Vasicek))

#___________________________Validation - risk neutrality_________________

deflators <-data.frame("deflator"=c(calcul$priceZCB[1],calcul$priceZCB[2], ESG_HW[2,1], ESG_HW[2,2],
  ESG_HW[2,3],ESG_HW[2,4], ESG_HW[3,5], ESG_HW[4,6], ESG_HW[3,7], ESG_HW[4,8], 
  ESG_HW[6,9], ESG_HW[11,10], ESG_HW[11,11], ESG_HW[11,12], ESG_HW[11,13]))
deflators$expected_price <-rep(deflators[1,1],15)
for (i in 2:15){
  deflators$expected_price[i] <- deflators[i,1]*deflators$expected_price[i-1]
}


install.package('vscDebugger')  
