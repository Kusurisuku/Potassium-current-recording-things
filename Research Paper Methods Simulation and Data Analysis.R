#simulate data

library(minpack.lm) #library for exponential fitting
library(ggplot2) #library for plotting
library(tidyr) #library for data format conversion from table data frame to a long list, which is easier for plotting
library(agricolae) #library for the Student–Newman–Keuls post hoc test

#Create a list of data to dissect, imaging this is a cell with ~100pF cellular capacitance
data_list <- lapply(1:200, function(i) {
t<- seq(0, 15000, by = 1) # time in ms
y_no_noise <- 1800*exp(-t/60) + 1500*exp(-t/800) + 1300*exp(-t/4000) + 400 # simulate a 3-exponential decay with errors
error_sd <- sqrt(y_no_noise) #poission error
y <- y_no_noise + rnorm(length(t), 0, sd=error_sd)
data.frame(t = t, y = y)
})

#create a list for results from each loop
result_2exp <- list()
result_3exp <- list()
result_4exp <- list()

nlc <- nls.control(maxiter = 1000) #increase the fitting attempts to 1000 from default

for (i in 1:length(data_list)) {
  
  y <- data_list[[i]]$y #y axis
  t <- data_list[[i]]$t #x axis
  
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

y_fit <- predict(fit_2exp)

# Plot other fitted components

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

# Generate fitted curve
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
result_2exp[[i]] <- list()
result_2exp[[i]]<- data.frame(
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

y_fit <- predict(fit_3exp)

# Plot other fitted components

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

# Generate fitted curve
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
result_3exp[[i]] <- list()
result_3exp[[i]]<- data.frame(
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

fit_4exp <- try(nlsLM(
  y ~ A1*exp(-t/tau1) + A2*exp(-t/tau2) + A3*exp(-t/tau3) + A4*exp(-t/tau4) + C,
  control=nlc,
  start = list(A1=1800, tau1=60, A2=1500, tau2=800, A3=700, tau3=4000, A4=1500, tau4=4000, C=40),
  lower = c(A1=0, tau1=0, A2=0, tau2=0, A3=0, tau3=0, A4=0, tau4=0, C=0),
  upper = c(A1=max_current, tau1=max_tau, A2=max_current, tau2=max_tau, A3=max_current, tau3=max_tau, A4=max_current, tau4=max_tau, C=max_current)
  ),
  silent=FALSE)
if (inherits(fit_4exp, "try-error")) next

y_fit <- predict(fit_4exp)

# Plot other fitted components

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

# Generate fitted curve
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
result_4exp[[i]] <- list()
result_4exp[[i]]<- data.frame(
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

result_4exp <- Filter(Negate(is.null), result_4exp) #remove Null results
result_2exp <- result_2exp[1:105] #correct length
result_3exp <- result_3exp[1:105] #correct length
result_4exp <- result_4exp[1:105] #correct length

#First we analyze the correlation coefficient R

r2 <- sapply(result_2exp, function(x) x$R2)
r3 <- sapply(result_3exp, function(x) x$R2)
r4 <- sapply(result_4exp, function(x) x$R2)
r_table <- data.frame("2exp"=r2, "3exp"=r3, "4exp"=r4, check.names = FALSE) #check.name remove annoying "X" in front of any titles start with a number

r_list <- pivot_longer(r_table, cols = everything(), names_to = "Exponential", values_to = "R")
ggplot(r_list, aes(x = Exponential, y = R)) +
  geom_jitter(width = 0.05, size = 1.5, alpha = 0.5) +
  theme_minimal() +
  labs(y = "Correlation Coefficient (R)", title = "Fitting Quality for Each Exponential Model") +
  theme(legend.position = "none")

summary(aov(R ~ Exponential, data = r_list))
print(SNK.test(aov(R ~ Exponential, data = r_list),"Exponential", group=TRUE))

#Then the Reduced chi-squared

rcs2 <- sapply(result_2exp, function(x) x$red_chi_2exp)
rcs3 <- sapply(result_3exp, function(x) x$red_chi_3exp)
rcs4 <- sapply(result_4exp, function(x) x$red_chi_4exp)
rcs_table <- data.frame("2exp"=rcs2, "3exp"=rcs3, "4exp"=rcs4, check.names = FALSE) #check.name into false remove annoying "X" in front of any titles start with a number

rcs_list <- pivot_longer(rcs_table, cols = everything(), names_to = "Exponential", values_to = "RCS")
ggplot(rcs_list, aes(x = Exponential, y = RCS)) +
  ylim(0, 3) +
  geom_jitter(width = 0.05, size = 1.5, alpha = 0.5) +
  theme_minimal() +
  labs(y = "Reduced Chi-squared X^2", title = "Fitting Quality for Each Exponential Model") +
  theme(legend.position = "none")

summary(aov(RCS ~ Exponential, data = rcs_list))
print(SNK.test(aov(RCS ~ Exponential, data = rcs_list),"Exponential", group=TRUE))

#F-distribution test is a good way to test the necessity of adding an extra component

SSE_2 <- as.numeric(sapply(result_2exp, function(x) x$SSE_2exp))
SSE_3 <- as.numeric(sapply(result_3exp, function(x) x$SSE_3exp))
SSE_4 <- as.numeric(sapply(result_4exp, function(x) x$SSE_4exp))
k_2 <- 5
k_3 <- 7
k_4 <- 9
n <- length(t)

T_2to3 <- ((SSE_2 - SSE_3) / SSE_3) * ((n - k_3) / 2)
T_3to4 <- ((SSE_3 - SSE_4) / SSE_4) * ((n - k_4) / 3)
T_table <- data.frame("2-3exp"=T_2to3, "3-4exp"=T_3to4, check.names = FALSE)
p_2to3 <- 1 - pf(T_2to3, df1 = k_3 - k_2, df2 = n - k_3)
p_3to4 <- 1 - pf(T_3to4, df1 = k_4 - k_3, df2 = n - k_4)
p_table <- data.frame("2-3exp"=p_2to3, "3-4exp"=p_3to4, check.names = FALSE)

T_list <- pivot_longer(T_table, cols = everything(), names_to = "Exponential", values_to = "TV")
p_list <- pivot_longer(p_table, cols = everything(), names_to = "Exponential", values_to = "PV")
ggplot(T_list, aes(x = Exponential, y = TV)) +
  geom_jitter(width = 0.05, size = 1.5, alpha = 0.5) +
  theme_minimal() +
  labs(y = "T Value", title = "Fitting Quality for Each Exponential Model") +
  theme(legend.position = "none")

t.test(TV ~ Exponential, data = T_list)
t.test(PV ~ Exponential, data = p_list)


#Finally, we plot the dissected Tau from each exponential fits

t4 <- c(sapply(result_4exp, function(x) x$Tau1), sapply(result_4exp, function(x) x$Tau2), sapply(result_4exp, function(x) x$Tau3), sapply(result_4exp, function(x) x$Tau4))
t3 <- c(sapply(result_4exp, function(x) x$Tau1), sapply(result_4exp, function(x) x$Tau2), sapply(result_4exp, function(x) x$Tau3), rep((-1), (length(t4)/4)))
t2 <- c(sapply(result_4exp, function(x) x$Tau1), sapply(result_4exp, function(x) x$Tau2), rep((-1), (length(t4)/4*2)))
t_table <- data.frame("2exp"=t2, "3exp"=t3, "4exp"=t4, check.names = FALSE) #check.name remove annoying "X" in front of any titles start with a number

t_list <- pivot_longer(t_table, cols = everything(), names_to = "Exponential", values_to = "Tau")
ggplot(t_list, aes(x = Exponential, y = Tau)) +
  scale_y_log10() +
  geom_jitter(width = 0.05, size = 1.5, alpha = 0.5) +
  theme_minimal() +
  labs(y = "Inactivation Time constant Tau (ms)", title = "Fitting Quality for Each Exponential Model") +
  theme(legend.position = "none")