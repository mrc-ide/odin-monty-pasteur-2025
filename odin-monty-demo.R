## Install packages needed (this will install odin2/dust2/monty, and additional
## packages used in this script)
install.packages(
  c("odin2", "dust2", "monty", "decor", "pkgload", "posterior", "brio"),
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

## Check if you have RTools on Windows,
## XCode command line tools on macOS or a 
## functioning C++ toolchain on Linux
pkgbuild::check_build_tools(debug = TRUE)

## Read in incidence data; see data/incidence.csv for original
data <- read.csv("data/incidence.csv")
plot(data, pch = 19, col = "red")

set.seed(1)

## Load packages

library(odin2)
library(dust2)
library(monty)


## Compiling the ODE model with `odin2`

sir_ode <- odin({
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
})



## Running the model with `dust2`

sys <- dust_system_create(sir_ode, pars = list())
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)

# The output has dimensions number of states x number of timepoints
dim(y)



## Unpacking states

y <- dust_unpack_state(sys, y)

plot(t, y$I, type = "l", xlab = "Time", ylab = "I")



## Compiling the stochastic SIR model

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
})



## Running a single simulation

sys <- dust_system_create(sir, pars = list(), dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

plot(t, y$I, type = "l", xlab = "Time", ylab = "Infected population")



## Running multiple simulations

sys <- dust_system_create(sir, pars = list(), n_particles = 50,
                          dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)
matplot(t, t(y$I), type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infected population")



## Calculating incidence with `zero_every`

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
})



## Incidence accumulates then resets

sys <- dust_system_create(sir, pars = list(), dt = 1 / 128)
dust_system_set_state_initial(sys)
t <- seq(0, 20, by = 1 / 128)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

plot(t[t %% 1 == 0], y$incidence[t %% 1 == 0], type = "l",
     ylab = "Infection incidence", xlab = "Time")
lines(t, y$incidence, col = "red")



## Time-varying inputs: using time

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  seed <- if (time == seed_time) seed_size else 0
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- min(seed + Binomial(S, p_SI), S)
  n_IR <- Binomial(I, p_IR)
  
  initial(S) <- N
  initial(I) <- 0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  seed_time <- parameter()
  seed_size <- parameter()
})


pars <- list(seed_time = 10, seed_size = 15)
sys <- dust_system_create(sir, pars = pars, dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

plot(t, y$I, type = "l", xlab = "Time", ylab = "Infected population")



## Time-varying inputs: using interpolate

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  beta <- interpolate(beta_time, beta_value, "constant")
  beta_time <- parameter()
  beta_value <- parameter()
  dim(beta_time, beta_value) <- parameter(rank = 1)
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  gamma <- parameter(0.1)
})



## Time-varying inputs: using `interpolate`

pars <- list(beta_time = c(0, 30), beta_value = c(0.2, 1))
sys <- dust_system_create(sir, pars = pars, dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

plot(t, y$I, type = "l", xlab = "Time", ylab = "Infected population")



## Using arrays
## A model with age stratification
sir_age <- odin({
  # Equations for transitions between compartments by age group
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  update(incidence) <- incidence + sum(n_SI)
  
  # Individual probabilities of transition:
  p_SI[] <- 1 - exp(-lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  # Calculate force of infection
  
  # age-structured contact matrix: m[i, j] is mean number of contacts an
  # individual in group i has with an individual in group j per time unit
  
  m <- parameter()
  
  # here s_ij[i, j] gives the mean number of contacts and individual in group
  # i will have with the currently infectious individuals of group j
  s_ij[, ] <- m[i, j] * I[j]
  
  # lambda[i] is the total force of infection on an individual in group i 
  lambda[] <- beta * sum(s_ij[i, ])
  
  # Draws from binomial distributions for numbers changing between
  # compartments:
  n_SI[] <- Binomial(S[i], p_SI[i])
  n_IR[] <- Binomial(I[i], p_IR)
  
  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  # User defined parameters - default in parentheses:
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  # Dimensions of arrays
  n_age <- parameter()
  dim(S, S0, n_SI, p_SI) <- n_age
  dim(I, I0, n_IR) <- n_age
  dim(R) <- n_age
  dim(m, s_ij) <- c(n_age, n_age)
  dim(lambda) <- n_age
})

pars <- list(S0 = c(990, 1000),
             I0 = c(10, 0),
             m = array(c(1.8, 0.4, 0.4, 1.2) / 2000, dim = c(2, 2)),
             beta = 0.2,
             gamma = 0.1,
             n_age = 2)

sys <- dust_system_create(sir_age, pars = pars, dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

matplot(t, t(y$I), type = "l", lty = 1, col = c("red", "blue"),
        xlab = "Time", ylab = "Infected population")
legend("topright", c("children", "adults"), col = c("red", "blue"), lty = 1)



## Arrays: model with age and vaccine

sir_age_vax <- odin({
  # Equations for transitions between compartments by age group
  update(S[, ]) <- new_S[i, j]
  update(I[, ]) <- I[i, j] + n_SI[i, j] - n_IR[i, j]
  update(R[, ]) <- R[i, j] + n_IR[i, j]
  update(incidence) <- incidence + sum(n_SI)
  
  # Individual probabilities of transition:
  p_SI[, ] <- 1 - exp(-rel_susceptibility[j] * lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  p_vax[, ] <- 1 - exp(-eta[i, j] * dt)
  
  # Force of infection
  m <- parameter() # age-structured contact matrix
  s_ij[, ] <- m[i, j] * sum(I[j, ])
  lambda[] <- beta * sum(s_ij[i, ])
  
  # Draws from binomial distributions for numbers changing between
  # compartments:
  n_SI[, ] <- Binomial(S[i, j], p_SI[i, j])
  n_IR[, ] <- Binomial(I[i, j], p_IR)
  
  # Nested binomial draw for vaccination in S
  # Assume you cannot move vaccine class and get infected in same step
  n_S_vax[, ] <- Binomial(S[i, j] - n_SI[i, j], p_vax[i, j])
  new_S[, 1] <- S[i, j] - n_SI[i, j] - n_S_vax[i, j] + n_S_vax[i, n_vax]
  new_S[, 2:n_vax] <- S[i, j] - n_SI[i, j] - n_S_vax[i, j] + n_S_vax[i, j - 1]
  
  initial(S[, ]) <- S0[i, j]
  initial(I[, ]) <- I0[i, j]
  initial(R[, ]) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  # User defined parameters - default in parentheses:
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.0165)
  gamma <- parameter(0.1)
  eta <- parameter()
  rel_susceptibility <- parameter()
  
  # Dimensions of arrays
  n_age <- parameter()
  n_vax <- parameter()
  dim(S, S0, n_SI, p_SI) <- c(n_age, n_vax)
  dim(I, I0, n_IR) <- c(n_age, n_vax)
  dim(R) <- c(n_age, n_vax)
  dim(m, s_ij) <- c(n_age, n_age)
  dim(lambda) <- n_age
  dim(eta) <- c(n_age, n_vax)
  dim(rel_susceptibility) <- c(n_vax)
  dim(p_vax, n_S_vax, new_S) <- c(n_age, n_vax)
})


set.seed(1)


## Model with likelihood

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  cases <- data()
  cases ~ Poisson(incidence)
})



## Calculating likelihood

filter <- dust_filter_create(sir, data = data, time_start = 0,
                             n_particles = 200, dt = 0.25)
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))



## Filtered trajectories

dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2),
                    save_trajectories = TRUE)
y <- dust_likelihood_last_trajectories(filter)
y <- dust_unpack_state(filter, y)
matplot(data$time, t(y$incidence), type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
points(data, pch = 19, col = "red")



## Parameter packers

packer <- monty_packer(c("beta", "gamma"))
packer

packer$unpack(c(0.2, 0.1))

packer$pack(c(beta = 0.2, gamma = 0.1))

packer <- monty_packer(c("beta", "gamma"), fixed = list(I0 = 5))
packer$unpack(c(0.2, 0.1))

packer <- monty_packer(array = c(beta = 3, gamma = 3))
packer
packer$unpack(c(0.2, 0.21, 0.22, 0.1, 0.11, 0.12))



## Priors

prior <- monty_dsl({
  beta ~ Exponential(mean = 0.5)
  gamma ~ Exponential(mean = 0.3)
})
prior

monty_model_density(prior, c(0.2, 0.1))

dexp(0.2, 1 / 0.5, log = TRUE) + dexp(0.1, 1 / 0.3, log = TRUE)



## From a dust filter to a monty model

filter

packer <- monty_packer(c("beta", "gamma"))
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE)
likelihood



## Posterior from likelihood and prior

posterior <- likelihood + prior
posterior



## Create a sampler

vcv <- diag(2) * 0.2
vcv

sampler <- monty_sampler_random_walk(vcv)
sampler



## Let's sample!

samples <- monty_sample(posterior, sampler, 1000, n_chains = 3)
samples



## The result: diagnostics

samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)



## The results: parameters

bayesplot::mcmc_scatter(samples_df)



## The result: traceplots

bayesplot::mcmc_trace(samples_df)



## The result: density over time

matplot(drop(samples$density), type = "l", lty = 1)
matplot(drop(samples$density[-(1:100), ]), type = "l", lty = 1)



## Better mixing

vcv <- matrix(c(0.01, 0.005, 0.005, 0.005), 2, 2)
sampler <- monty_sampler_random_walk(vcv)
samples <- monty_sample(posterior, sampler, 2000, initial = samples,
                        n_chains = 4)
matplot(samples$density, type = "l", lty = 1)



## Better mixing: the results

samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)

bayesplot::mcmc_scatter(samples_df)

bayesplot::mcmc_trace(samples_df)



## Trajectories

trajectories <- dust_unpack_state(filter,
                                  samples$observations$trajectories)
matplot(data$time, trajectories$incidence[, , 1], type = "l", lty = 1,
        col = "#00000044", xlab = "Time", ylab = "Infection incidence")
points(data, pch = 19, col = "red")

dim(samples$observations$trajectories)



## Saving a subset of trajectories

likelihood <- dust_likelihood_monty(filter, packer, 
                                    save_trajectories = c("I", "incidence"))
posterior <- likelihood + prior
samples2 <- monty_sample(posterior, sampler, 100, initial = samples)
dim(samples2$observations$trajectories)



## Thinning

samples <- monty_samples_thin(samples,
                              burnin = 500,
                              thinning_factor = 2)


set.seed(1)


## Fitting in deterministic mode

unfilter <- dust_unfilter_create(sir, data = data, time_start = 0, dt = 0.25)

dust_likelihood_run(unfilter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(unfilter, list(beta = 0.4, gamma = 0.2))

likelihood <- dust_likelihood_monty(unfilter, packer, save_trajectories = TRUE)
posterior <- likelihood + prior
samples_det <- monty_sample(posterior, sampler, 1000, n_chains = 4)
samples_det <- monty_samples_thin(samples_det,
                                  burnin = 500,
                                  thinning_factor = 2)



## Stochastic v deterministic comparison

y <- dust2::dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(20, 1000))
matplot(data$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75))
points(data, pch = 19, col = "red")

y <- dust2::dust_unpack_state(filter, samples_det$observations$trajectories)
incidence <- array(y$incidence, c(20, 1000))
matplot(data$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75))
points(data, pch = 19, col = "red")

pars_stochastic <- array(samples$pars, c(2, 500))
pars_deterministic <- array(samples_det$pars, c(2, 500))
plot(pars_stochastic[1, ], pars_stochastic[2, ], ylab = "gamma", xlab = "beta",
     pch = 19, col = "blue")
points(pars_deterministic[1, ], pars_deterministic[2, ], pch = 19, col = "red")
legend("bottomright", c("stochastic fit", "deterministic fit"), pch = c(19, 19), 
       col = c("blue", "red"))



## Projections and counterfactuals

data <- read.csv("data/schools.csv")
plot(data, pch = 19, col = "red")

sis <- odin({
  update(S) <- S - n_SI + n_IS
  update(I) <- I + n_SI - n_IS
  update(incidence) <- incidence + n_SI
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(incidence, zero_every = 1) <- 0
  
  schools <- interpolate(schools_time, schools_open, "constant")
  schools_time <- parameter()
  schools_open <- parameter()
  dim(schools_time, schools_open) <- parameter(rank = 1)
  
  beta <- ((1 - schools) * (1 - schools_modifier) + schools) * beta0
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IS <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IS <- Binomial(I, p_IS)
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta0 <- parameter(0.2)
  gamma <- parameter(0.1)
  schools_modifier <- parameter(0.6)
  
  cases <- data()
  cases ~ Poisson(incidence)
})

schools_time <- c(0, 50, 60, 120, 130, 170, 180)
schools_open <- c(1,  0,  1,   0,   1,   0,   1)


## Fitting to the SIS model

packer <- monty_packer(c("beta0", "gamma", "schools_modifier"),
                       fixed = list(schools_time = schools_time,
                                    schools_open = schools_open))

filter <- dust_filter_create(sis, time_start = 0, dt = 1,
                             data = data, n_particles = 200)

prior <- monty_dsl({
  beta0 ~ Exponential(mean = 0.3)
  gamma ~ Exponential(mean = 0.1)
  schools_modifier ~ Uniform(0, 1)
})

vcv <- diag(c(2e-4, 2e-4, 4e-4))
sampler <- monty_sampler_random_walk(vcv)


likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE,
                                    save_state = TRUE, save_snapshots = 60)

posterior <- likelihood + prior

samples <- monty_sample(posterior, sampler, 500, initial = c(0.3, 0.1, 0.5),
                        n_chains = 4)
samples <- monty_samples_thin(samples, burnin = 100, thinning_factor = 8)



## Fit to data

y <- dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(150, 200))
matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
points(data, pch = 19, col = "red")



## Running projection using the end state

state <- array(samples$observations$state, c(3, 200))
pars <- array(samples$pars, c(3, 200))
pars <- lapply(seq_len(200), function(i) packer$unpack(pars[, i]))

sys <- dust_system_create(sis, pars, n_groups = length(pars), dt = 1)

dust_system_set_state(sys, state)
t <- seq(150, 200)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence", xlim = c(0, 200))
matlines(t, t(y$incidence), col = "blue")
points(data, pch = 19, col = "red")



## Running counterfactual using the snapshot

snapshot <- array(samples$observations$snapshots, c(3, 200))
pars <- array(samples$pars, c(3, 200))
f <- function(i) {
  p <- packer$unpack(pars[, i])
  p$schools_time <- c(0, 50, 130, 170, 180)
  p$schools_open <- c(1, 0, 1, 0, 1)
  p
}
pars <- lapply(seq_len(200), f)
sys <- dust_system_create(sis, pars, n_groups = length(pars), dt = 1)

dust_system_set_state(sys, snapshot)
t <- seq(60, 150)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
matlines(t, t(y$incidence), col = "blue")
points(data, pch = 19, col = "red")


set.seed(1)


## Advanced `monty`: parallel tempering
## Mixture of normals model
ex_mixture <- function(x) log(dnorm(x, mean = 5) / 2 + dnorm(x, mean = -5) / 2)
likelihood <- monty_model_function(ex_mixture, allow_multiple_parameters = TRUE)

## Prior with wider variance
prior <- monty_dsl({
  x ~ Normal(0, 10)
})

## posterior is the product of the two (the sum on a log-scale)
posterior <- likelihood + prior
x <- seq(-10, 10, length.out = 1001)
y <- exp(posterior$density(rbind(x)))
plot(x, y / sum(y) / diff(x)[[1]], col = "red", type = 'l', ylab = "density")



## Try random walk Metropolis-Hastings sampler
vcv <- matrix(1.5)
sampler_rw <- monty_sampler_random_walk(vcv = vcv)
res_rw <- monty_sample(posterior, sampler_rw, n_steps = 1000)

## Only one mode has been explored
plot(res_rw$pars[1, , ], ylim = c(-8, 8), ylab = "Value")
abline(h = c(-5, 5), lty = 3, col = "red")

## Multiple chains does not solve the issue

res_rw3 <- monty_sample(posterior, sampler_rw, n_steps = 1000, n_chains = 3,
                        initial = rbind(c(-5, 0, 5)))
matplot(res_rw3$pars[1, , ], type = "p", pch = 21, ylim = c(-8, 8),
        ylab = "Value", col = c("blue", "green", "purple"))
abline(h = c(-5, 5), lty = 3, col = "red")



## Use our RW sampler within the parallel tempering sampler
vcv <- matrix(1.5)
sampler_rw <- monty_sampler_random_walk(vcv = vcv)
s <- monty_sampler_parallel_tempering(sampler_rw, n_rungs = 10)
res_pt <- monty_sample(posterior, s, 1000, n_chains = 1)

plot(res_pt$pars[1, , ], ylim = c(-8, 8), ylab = "Value")
abline(h = c(-5, 5), lty = 3, col = "red")
