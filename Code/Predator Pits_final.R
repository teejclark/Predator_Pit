# Code to simulate stochastic density dependent population growth of prey with predation

#### T.J. Clark et al. - Stochastic predator pits...

# load packages
library(MASS)
library(logitnorm)

# Model for prey population growth:
# N(t+1)/N(t) = Exp[r_max*(1-(N(t)/K)^theta)+E(t)]-InvLogit{Logit[P(t)]+F(t)}*A*(1+R)  ; E(t)~N(0,sigmaE; F(t)~N(0,sigmaP))

K <- 20 # carrying capacity (#elk/km^2); "Low" = 5; "High" = 20

r_max <- 0.28 # r=max of 0.28 equals lambda_max of 1.32; 0.97 for r-selected
# r_max = 0.97  #r=max of 0.97 equals 2.64 for r-selected

theta <- 4 # shape parameter for density dependent response in theta-logistic model (4 for elk; 0.3 for r-selected)
# theta = 0.3   #shape parameter for density dependent response in theta-logistic model (4 for elk; 0.3 for r-selected)

A <- 0.7 # proportion of predation that is additive; norm=0.7; cancel=1
R <- exp(r_max) - 1 # elk reproductive rate (per capita recruitement rate max. growth rate)
sigmaE <- 0.1 # amount of environmental stochasticity in prey population growth rate  (0.1)
sigmaP <- 1.75 #amount of variation in stochastic predation rate: "Low" = 1.1; "High" = 2.3
Sigma <- matrix(c(1, .7, .7, 1), 2, 2) # correlation between sigmaE and sigmaP

# generate stochasticity

# Model for predation rate (proportion of prey killed/year):
# P(t) = Psi(t)*W(t)/[N(t)*1000]

# Model for functional response(kill rate = Psi(t); #prey killed/predator/year)
# Psi(t)=[alpha_0*N(t)^(alpha_2+1)]/[alpha_1+N(t)^(alpha_2+1)]

alpha_0 <- 20 # maximum kill rate; 20
alpha_1 <- 1 # rate of increase in kill rate  ; 0.4
alpha_2 <- 0 # shape parameter for functional response; 0 = Type II

# Model for numerical response (W(t) = prey density in #/1000km^2)
# W(t) = delta_0+[(delta_1-delta_0)*N(t)]/[delta_2+N(t)]

delta_0 <- 5 # minimum predator density at prey population; 5
delta_1 <- 35 # maximum predator density; 35
delta_2 <- 1 # rate of increase in predator density; 1

InitAbund <- 3 # orig val = 3
NumReps <- 1000
NumTimeSteps <- 500

# generate stochasticity
temp <- mvrnorm(n = NumReps * NumTimeSteps, mu = c(0, 0), Sigma)

# End SetUp

#-----------------------------------------------------------------------------

# Create empty matrices
N_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps + 1))
colnames(N_t) <- paste(rep("t_", NumTimeSteps + 1), 0:NumTimeSteps, sep = "")
W_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(W_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
Psi_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(Psi_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
P_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(P_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
Pred_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(Pred_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
r_t <- as.data.frame(matrix(-Inf, NumReps, NumTimeSteps))
colnames(r_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
lambda_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(lambda_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")
Stoch_P_t <- as.data.frame(matrix(0, NumReps, NumTimeSteps))
colnames(lambda_t) <- paste(rep("t_", NumTimeSteps), 0:(NumTimeSteps - 1), sep = "")

N_t[, 1] <- InitAbund


# matrix to save extinction values
#mats <- matrix(numeric(length(sigmaP) * length(alpha_0)), nrow = length(sigmaP))
mats <- array(NA,
  dim = c(NumReps, NumTimeSteps, length(sigmaP)),
  dimnames = list(
    c(paste(rep("Rep_", NumReps), 1:NumReps, sep = "")),
    c(paste(rep("t_", NumTimeSteps), 1:NumTimeSteps, sep = "")),
    c(1:length(sigmaP))
  )
)


# start the simulation
# for (i in 1:length(sigmaP)){
#  for (j in 1:length(alpha_0)){

#for (i in 1:length(sigmaP)){
for (t in 1:NumTimeSteps) {

  # stochastic elements
  # f_t = matrix(0,NumReps,NumTimeSteps)
  # e_t = matrix(0,NumReps,NumTimeSteps)
  f_t <- matrix(temp[, 2] * sigmaP, NumReps, NumTimeSteps)
  e_t <- matrix(temp[, 1] * sigmaE, NumReps, NumTimeSteps)

  W_t[, t] <- (delta_0 + ((delta_1 - delta_0) * N_t[, t-1]) / (delta_2 + N_t[, t-1])) # numerical response
  Psi_t[, t] <- (alpha_0 * N_t[, t]^(alpha_2 + 1)) / (alpha_1 + N_t[, t]^(alpha_2 + 1)) # functional response
  P_t[, t] <- (Psi_t[, t] * W_t[, t] / (N_t[, t] * 1000)) # total response
  Pred_t[, t] <- invlogit(logit(P_t[, t] * A)+ f_t[, t]) # predation, including the additive effect

  r_t[, t] <- r_max * (1 - (N_t[, t] / K)^theta) + e_t[, t]
  r_t[which(r_t[,t]>r_t_max),t] = r_t_max

  lambda_t[, t] <- exp(r_t[, t])*(1-Pred_t[,t])
  lambda_t[which(lambda_t[, t] == 0), t] <- 0.01

  N_t[, t + 1] <- N_t[, t] * lambda_t[, t] # prey growth

  # limit how low the population can go
  #N_t[which(N_t[,t] < 0.5) , t] <- 0.5

  #mats[,t,i] = N_t[,t]


  # extinction
  #mats[i,j] = length(which(N_t[,t]<(1)))/NumReps
}

plot(density(N_t$t_500),
     main = "", xlim = c(0,K), lwd = 4, axes = T,
     xlab = "Abundance", ylab = "Density",cex.axis=2,cex.lab=2
)
abline(v=5,lwd=2,lty=2)

#### Calculate Properties of the Simulations

#### Stochastic sims...fractions of simulations ending at high equilibrium
# note - K = 12 is an arbitrary-ish number
length(which(N_t[,500] > 15))/NumReps # 15


# proportion near K
PropNearK <- length(which(N_t[, 40:NumTimeSteps] > (0.90 * K))) / (length(40:NumTimeSteps) * NumReps)
PropNearK

# Proportion in Pit
PropInPit <- length(which(N_t[, NumTimeSteps] < (0.2 * K))) / NumReps
PropInPit

# Proportion Extinct
PropExtinct <- length(which(N_t[, NumTimeSteps] < (0.5))) / NumReps
PropExtinct


### Graph Extinction Rates

# Matrix Plot
library(ggplot2)
library(reshape)
finalmat <- melt(mats)

windows()
ggplot(finalmat) +
  aes(x = X1, y = X2, z = value, fill = value) +
  geom_tile() +
  coord_equal() +
  # geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  labs(
    title = "",
    x = "Predation Stochasticity",
    y = "Max Predation Rate (elk/wolf/year)",
    z = "Extinction Likelihood (%)"
  ) +
  theme_bw()


################################################################################################
