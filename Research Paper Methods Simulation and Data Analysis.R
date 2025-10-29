#simulate data

#Create a list of data to dissect, imaging this is a cell with ~100pF cellular capacitance
data_list <- lapply(1:3, function(i) {
t<- seq(0, 15000, by = 1) # time in ms
y_no_noise <- 1800*exp(-t/60) + 1500*exp(-t/800) + 1300*exp(-t/4000) + 400 # simulate a 3-exponential decay with errors
error_sd <- sqrt(y_no_noise) #poission error
y <- y_no_noise + rnorm(length(t), 0, sd=error_sd)
data.frame(t = t, y = y)
})

max_current <- max(y) #for determine current upper boundary
max_tau <- max(t)# for determine tau upper boundary

library(minpack.lm)

for (i in 1:length(data_list)) {
  
  data <- data_list[[i]]

# 2 exponential fitting the simulated data with initial guesses and boundaries
fit_2exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + C,
  data = data,
  start = list(A1=4000, tau1=50, A2=1500, tau2=500, C=40),
  lower = c(A1=0, tau1=0, A2=0, tau2=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_2exp)

# Plot other fitted components

res <- resid(fit_2exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
R2

summary(fit_2exp)
coef(fit_2exp)

params <- coef(fit_2exp)
stats <- summary(fit_2exp)
dof <- stats$df[2]
chi_2exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_2exp <- chi_2exp/dof
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


# 3 exponential fitting the simulated data with initial guesses and boundaries
fit_3exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + A3*exp(-t/tau3) + C,
  data = data,
  start = list(A1=4000, tau1=50, A2=1500, tau2=500, A3=500, tau3=4000, C=40),
  lower = c(A1=0, tau1=0, A2=0, tau2=0, A3=0, tau3=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, A3=max_current, tau3=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_3exp)

# Plot other fitted components

res <- resid(fit_3exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
R2

summary(fit_3exp)
coef(fit_3exp)

params <- coef(fit_3exp)
stats <- summary(fit_3exp)
dof <- stats$df[2]
chi_3exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_3exp <- chi_3exp/dof
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

# 4 exponential fitting the simulated data with initial guesses and boundaries
fit_4exp <- nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + A3*exp(-t/tau3) + A4*exp(-t/tau4) + C,
  data = data,
  start = list(A1=2000, tau1=50, A2=1500, tau2=500, A3=1500, tau3=2000, A4=1500, tau4=5000, C=40),
  lower = c(A1=0, tau1=0, A2=0, tau2=0, A3=0, tau3=0, A4=0, tau4=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, A3=max_current, tau3=max_tau, A4=max_current, tau4=max_tau, C=max_current)
)


# Generate fitted curve
y_fit <- predict(fit_4exp)

# Plot other fitted components

res <- resid(fit_4exp)
RSS <- sum(res^2)
TSS <- sum((y - mean(y))^2)
R2 <- 1 - RSS/TSS
R2

summary(fit_4exp)
coef(fit_4exp)

params <- coef(fit_4exp)
stats <- summary(fit_4exp)
dof <- stats$df[2]
chi_4exp <- sum(((y-y_fit)^2)/(y_fit))
red_chi_4exp <- chi_4exp/dof
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
}