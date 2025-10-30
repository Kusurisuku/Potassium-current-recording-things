#simulate data

#Create a list of data to dissect, imaging this is a cell with ~100pF cellular capacitance
data_list <- lapply(1:3, function(i) {
t<- seq(0, 15000, by = 1) # time in ms
y_no_noise <- 1800*exp(-t/60) + 1500*exp(-t/800) + 1300*exp(-t/4000) + 400 # simulate a 3-exponential decay with errors
error_sd <- sqrt(y_no_noise) #poission error
y <- y_no_noise + rnorm(length(t), 0, sd=error_sd)
data.frame(t = t, y = y)
})


library(minpack.lm)
nlc <- nls.control(maxiter = 1000)

#create a list for results from each loop
result_list <- list()
result_list[["2exp"]] <- list()
result_list[["3exp"]] <- list()
result_list[["4exp"]] <- list()


for (i in 1:length(data_list)) {
  
  y <- data_list[[i]]$y
  t <- data_list[[i]]$t
  
  max_current <- max(y) #for determine current upper boundary
  max_tau <- max(t)# for determine tau upper boundary

# 2 exponential fitting the simulated data with initial guesses and boundaries
fit_2exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + C,
  start = list(A1=4000, tau1=50, A2=1500, tau2=500, C=40),
  control=nlc,
  lower = c(A1=0, tau1=0, A2=0, tau2=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_2exp)

# Plot other fitted components

summary(fit_2exp)
coef(fit_2exp)

res <- resid(fit_2exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
params <- coef(fit_2exp)
stats <- summary(fit_2exp)
dof <- stats$df[2]
chi_2exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_2exp <- chi_2exp/dof
SSE_2exp = sum(((y-y_fit)^2))
A1 <- params["A1"]; tau1 <- params["tau1"]
A2 <- params["A2"]; tau2 <- params["tau2"]; C <- params["C"]

comp1 <- A1*exp(-t/tau1)
comp2 <- A2*exp(-t/tau2)
steady <- rep(params["C"],length(t))

plot(t, y, type="p", main="2-exponential component analysis", xlab="Time (ms)", ylab="Current (pA)",ylim = c(-5, max_current+100))
lines(t, comp1, col="blue", lwd=2)
lines(t, comp2, col="green", lwd=2)
lines(t, y_fit, col="red", lwd=2)
lines(t, res, col="orange", lwd=0.5)
lines(t, steady, col="yellow", lwd=2)
legend("topright", legend=c("Data","Total Fit","τ1","τ2","Residule","Steady-state current"),
       col=c("black","red","blue","green","orange","yellow"), lty=1)

#write parameters from 2-exp fit into the result list
result_list[["2exp"]][[i]] <- list()
result_list[['2exp']][[i]]<- data.frame(
  ID = paste0("Fit_", i,"-2_exp"),
  Amp1 = A1,
  Amp2 = A2,
  Tau1 = tau1,
  Tau2 = tau2,
  SSE_2exp = SSE_2exp,
  red_chi_2exp = red_chi_2exp,
  R2 = R2
)

# 3 exponential fitting the simulated data with initial guesses and boundaries
fit_3exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + A3*exp(-t/tau3) + C,
  start = list(A1=4000, tau1=50, A2=1500, tau2=500, A3=500, tau3=4000, C=40),
  control=nlc,
  lower = c(A1=0, tau1=0, A2=0, tau2=0, A3=0, tau3=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, A3=max_current, tau3=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_3exp)

# Plot other fitted components

summary(fit_3exp)
coef(fit_3exp)

res <- resid(fit_3exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
params <- coef(fit_3exp)
stats <- summary(fit_3exp)
dof <- stats$df[2]
chi_3exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_3exp <- chi_3exp/dof
SSE_3exp = sum(((y-y_fit)^2))
A1 <- params["A1"]; tau1 <- params["tau1"]
A2 <- params["A2"]; tau2 <- params["tau2"]
A3 <- params["A3"]; tau3 <- params["tau3"]; C <- params["C"]

comp1 <- A1*exp(-t/tau1)
comp2 <- A2*exp(-t/tau2)
comp3 <- A3*exp(-t/tau3)
steady <- rep(params["C"],length(t))

plot(t, y, type="p", main="3-exponential component analysis", xlab="Time (ms)", ylab="Current (pA)",ylim = c(-5, max_current+100))
lines(t, comp1, col="blue", lwd=2)
lines(t, comp2, col="green", lwd=2)
lines(t, comp3, col="purple", lwd=2)
lines(t, y_fit, col="red", lwd=2)
lines(t, res, col="orange", lwd=0.5)
lines(t, steady, col="yellow", lwd=2)
legend("topright", legend=c("Data","Total Fit","τ1","τ2","τ3","Residule","Steady-state current"),
       col=c("black","red","blue","green","purple","orange","yellow"), lty=1)

#write parameters from 3-exp fit into the result list
result_list[["3exp"]][[i]] <- list()
result_list[["3exp"]][[i]]<- data.frame(
  ID = paste0("Fit_", i,"-3_exp"),
  Amp1 = A1,
  Amp2 = A2,
  Amp3 = A3,
  Tau1 = tau1,
  Tau2 = tau2,
  Tau3 = tau3,
  SSE_3exp = SSE_3exp,
  red_chi_3exp = red_chi_3exp,
  R2 = R2
)

# 4 exponential fitting the simulated data with initial guesses and boundaries
fit_4exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + A3*exp(-t/tau3) + A4*exp(-t/tau4) + C,
  control=nlc,
  start = list(A1=2000, tau1=50, A2=1500, tau2=600, A3=1000, tau3=2000, A4=1300, tau4=3000, C=40),
  lower = c(A1=0, tau1=0, A2=0, tau2=0, A3=0, tau3=0, A4=0, tau4=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, A3=max_current, tau3=max_tau, A4=max_current, tau4=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_4exp)

# Plot other fitted components

summary(fit_4exp)
coef(fit_4exp)

res <- resid(fit_4exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
params <- coef(fit_4exp)
stats <- summary(fit_4exp)
dof <- stats$df[2]
chi_4exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_4exp <- chi_4exp/dof
SSE_4exp = sum(((y-y_fit)^2))
A1 <- params["A1"]; tau1 <- params["tau1"]
A2 <- params["A2"]; tau2 <- params["tau2"]
A3 <- params["A2"]; tau3 <- params["tau3"]
A4 <- params["A4"]; tau4 <- params["tau4"]
C <- params["C"]

comp1 <- A1*exp(-t/tau1)
comp2 <- A2*exp(-t/tau2)
comp3 <- A3*exp(-t/tau3)
comp4 <- A4*exp(-t/tau4)
steady <- rep(params["C"],length(t))

plot(t, y, type="p", main="4-exponential component analysis", xlab="Time (ms)", ylab="Current (pA)",ylim = c(-5, max_current+100))
lines(t, comp1, col="blue", lwd=2)
lines(t, comp2, col="green", lwd=2)
lines(t, comp3, col="purple", lwd=2)
lines(t, comp4, col="cyan", lwd=2)
lines(t, y_fit, col="red", lwd=2)
lines(t, res, col="orange", lwd=0.5)
lines(t, steady, col="yellow", lwd=2)
legend("topright", legend=c("Data","Total Fit","τ1","τ2","τ3","τ4", "Residule","Steady-state current"),
       col=c("black","red","blue","green","purple","cyan","orange","yellow"), lty=1)

#write parameters from 4-exp fit into the result list
result_list[["4exp"]][[i]] <- list()
result_list[["4exp"]][[i]]<- data.frame(
  ID = paste0("Fit_", i,"-4_exp"),
  Amp1 = A1,
  Amp2 = A2,
  Amp3 = A3,
  Amp4 = A4,
  Tau1 = tau1,
  Tau2 = tau2,
  Tau3 = tau3,
  Tau4 = tau4,
  SSE_4exp = SSE_4exp,
  red_chi_4exp = red_chi_4exp,
  R2 = R2
)
}

#analyze listed results