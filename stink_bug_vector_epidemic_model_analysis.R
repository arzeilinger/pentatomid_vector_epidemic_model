#### Compartmental epidemic models for Pentatomid vector-borne pathogen systems
#### Developed by Adam Zeilinger for:  
#### Chapter 13. Pentatomids as vectors of plant pathogens. In McPherson, J. (ed.), Biology of Invasive Stink Bugs and Related Species (Pentatomoidea) (201)

# Load packages
library(lattice)
library(deSolve)
library(data.table)
library(tidyr)


######################################################################################################
#### Symbiotic pathogen model
######################################################################################################

#### Numerical simulation of symbiotic model
# Equations for model
SymModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    # Defended patch
    dS <- a*I - (beta*S*V)/(I + S)
    dI <- (beta*S*V)/(I + S) - a*I
    dU <- bu*((K - U)/K)*U + mu*V - (alpha*I*U)/(I + S) - du*U
    dV <- bv*((K - V)/K)*V + (alpha*I*U)/(I + S) - mu*V - dv*V
    return(list(c(dS, dI, dU, dV)))
  })
}

## Function to run simulations and output infective host (I) and infectious vector (V) densities
SimFunc1 <- function(x){
  x <- as.numeric(x)
  Pars <- c(alpha = x[1], beta = x[2], 
            mu = x[3], a = x[4],
            bu = x[5], bv = x[6],
            du = x[7], dv = x[8],
            K = x[9])
  State <- c(S = S, I = I, U = U, V = V)
  Time <- seq(0, 2000, by = 2) # time steps to 2000 by 2
  model.out <- as.data.frame(ode(func = SymModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  model.dat <- data.frame(model.out[nrow(model.out),c("time", "S", "I", "U", "V")])
  return(c(x[3], x[6], model.dat))
}

#### Define non-varying parameter values
alpha <- 0.4
beta <- 0.4
a <- 0.1
bu <- 0.3
du <- 0.2
dv <- 0.2
K <- 250

#### State variable starting values
S = 100; E = 0; I = 0
U = 199; V = 1
N = S + I
M = U + V


#### Run multiple simulations and return only equilibria
#### Varying bv (transovarial transmission)
# Define varying parameters
mu <- 0
bv <- seq(0,1,by=0.025)
# Gather all parameter values
params <- cbind(alpha, beta, mu, a, bu, bv, du, dv, K)
# Run simulations
SymRun <- as.data.frame(rbindlist(apply(params, 1, SimFunc1)))
# Construct data set with varying bv
bvPerc <- data.frame(param = "bv",
                     var = SymRun[,2],
                     PercI = (SymRun$I/N)*100,
                     PercV = (SymRun$V/(SymRun$U+SymRun$V))*100)


#### Varying mu (recovery rate)
# Define varying parameters
mu <- seq(0,1,by=0.025)
bv <- 0
# Gather all parameter values
params <- cbind(alpha, beta, mu, a, bu, bv, du, dv, K)
# Run simulations
SymRun <- as.data.frame(rbindlist(apply(params, 1, SimFunc1)))
# Construct data set with varying mu -- vector recovery/persistence
MuPerc <- data.frame(param = "mu",
                     var = SymRun[,1],
                     PercI = (SymRun$I/N)*100,
                     PercV = (SymRun$V/(SymRun$U+SymRun$V))*100)


# Combine data sets and prepare for plotting
SymResults <- rbind(MuPerc,bvPerc)
PlotResults <- gather(SymResults, "State", "Value", 3:4)

# Plotting Figure 13.2
jpeg(file = "Symbiotic_model_plot.jpg",
     width = 76, height = 76, units = "mm", res = 600)
  plot.new()
  xyplot(Value ~ var|param, groups = State, data = PlotResults,
         scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 0.5, relation = "free",
                       x = list(limits = c(0,1), at = seq(0,1,0.5),
                                labels = list(seq(0,1,0.5),rep("",3))),
                       y = list(limits = c(0, 110), at = seq(0,100,50),
                                labels = list(seq(0,100,50), seq(0,100,50)))
         ),
         xlab = list("Parameter value", cex = 0.7), 
         ylab = list("Percent infected", cex = 0.7),
         layout = c(1,2), pch = c(16), cex = 1,
         type = 'l', aspect = 1, strip = FALSE,
         col = "black", lty = c(1,2), lwd = 2)
  mtext(c("A", "B"), side = 3, cex = 0.7, adj = rep(0.3, 2), padj = c(-4.7, 10))
dev.off()



################################################################################
#### Facultative pathogens model
################################################################################
#### Numerical simulation of Facultative model
# Equations for model
FacModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    # Defended patch
    dS <- a*I + a*E - beta*S
    dE <- beta*S - (delta*E*I)/N - a*E
    dI <- (delta*E*I)/N - a*I
    return(list(c(dS, dE, dI)))
  })
}

## Function to run simulations and output infective host (I) and infectious vector (V) densities
FacFunc <- function(x){
  x <- as.numeric(x)
  Pars <- c(beta = x[1], delta = x[2], a = x[3])
  State <- c(S = S, E = E, I = I)
  Time <- seq(0, 2000, by = 2) # time steps to 1500 by 2
  model.out <- as.data.frame(ode(func = FacModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  model.dat <- data.frame(model.out[nrow(model.out),c("time", "S", "E", "I")])
  return(c(x[2],model.dat))
}

#### Define parameter values
beta <- 0.4
delta <- seq(0,1,by=0.025)
a <- 0.1

params <- cbind(beta, delta, a)

#### State variable starting values
S = 99; E = 0; I = 1
N = S + E + I

#### Run multiple simulations and return only equilibria
FacRun <- as.data.frame(rbindlist(apply(params, 1, FacFunc)))

# Construct data set with varying beta
BetaPerc <- data.frame(param = "beta",
                        var = FacRun[,1],
                        PercI = (FacRun$I/N)*100)
# Construct data set with varying delta
DeltaPerc <- data.frame(param = "delta",
                        var = FacRun[,1],
                        PercI = (FacRun$I/N)*100)

# Combine data sets
FacResults <- rbind(BetaPerc,DeltaPerc)


# Plotting Figure 13.3
jpeg("Facultative_model_plot.jpg",
     width = 76, height = 76, units = "mm", res = 600)
  plot.new()
  xyplot(PercI ~ var|param, data = FacResults,
         scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 0.5, relation = "free",
                       x = list(limits = c(0,1), at = seq(0,1,0.5),
                                labels = list(seq(0,1,0.5),rep("",3))),
                       y = list(limits = c(0, 110), at = seq(0,100,50),
                                labels = list(seq(0,100,50), seq(0,100,50)))
         ),
         xlab = list("Parameter value", cex = 0.7), 
         ylab = list("Percent infected", cex = 0.7),
         layout = c(1,2), pch = c(16), cex = 1,
         type = 'l', aspect = 1, strip = FALSE,
         col = "black", lty = c(1,2), lwd = 2)
  mtext(c("A", "B"), side = 3, cex = 0.7, adj = rep(0.3, 2), padj = c(-4.7, 10))
dev.off()