# 1. Install and load necessary packages
# install.packages("ergm")
# install.packages("network")
library(ergm)
library(network)
library(parallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 2. Read the data
# Note: Adjust 'header' and 'row.names' depending on the exact structure of your CSVs.
# Assuming standard adjacency matrix format and a header row for demographics.
f1_matrix <- as.matrix(read.csv("Glasgow/f1.csv", header = FALSE)) 
demographics <- read.csv("Glasgow/demographic.csv", header = TRUE)

# 3. Prepare the network object
# The description says "Pupils were asked to name...", indicating a directed network.
glasgow_net1 <- network(f1_matrix, directed = TRUE)

# Add the gender attribute to the network nodes
# Gender is coded as 1=boy, 2=girl
set.vertex.attribute(glasgow_net1, "gender", demographics$gender)

# 4. Estimate the ERGM
# We include 'edges' (density) and 'nodematch' for gender homophily
model_1 <- ergm(glasgow_net1 ~ edges + nodematch("gender"))

# View the results
summary(model_1)

# 5. Compute the conditional probability
# Extract the estimated coefficients (theta values)
coefs <- coef(model_1)
theta_edges <- coefs["edges"]
theta_gender_match <- coefs["nodematch.gender"]

# Calculate the log-odds of a tie between two students of the SAME gender
# Change statistic for edges is 1, change statistic for matching gender is 1
log_odds_same_gender <- theta_edges + theta_gender_match

# Convert log-odds to probability
prob_same_gender <- exp(log_odds_same_gender) / (1 + exp(log_odds_same_gender))

# Print the probability
cat("Conditional probability of a tie between same-gender students:", prob_same_gender, "\n")

# 1. Read the additional control data
logdistance_matrix <- as.matrix(read.csv("Glasgow/logdistance.csv", header = FALSE))
alcohol <- read.csv("Glasgow/alcohol.csv", header = TRUE) 

# 2. Add the new attributes to your existing network object (glasgow_net1)
set.vertex.attribute(glasgow_net1, "age", demographics$age)

# Assuming alcohol.csv has a column named 'alcohol' (adjust if named differently)
set.vertex.attribute(glasgow_net1, "alcohol", alcohol$Wave.1) 

# 3. Specify and estimate the comprehensive ERGM
# We fix the decay parameters to 0.6 as suggested by your assignment note
ergm_model_2_params <- glasgow_net1 ~ edges + 
  # Endogenous Hypotheses
  mutual + 
  gwesp(0.2, fixed = TRUE) + 
  
  # Exogenous Hypotheses: Gender
  nodeofactor("gender", levels = 1) + 
  nodeifactor("gender", levels = 1) + 
  
  # Exogenous Hypotheses: Alcohol
  nodeocov("alcohol") + 
  nodeicov("alcohol") + 
  absdiff("alcohol") +
  
  # Basic Exogenous Controls
  edgecov(logdistance_matrix) + 
  absdiff("age") +

  # Basic Endogenous Controls (Stabilizers)
  gwidegree(0.6, fixed = TRUE) + 
  gwodegree(0.6, fixed = TRUE)

noisy_model_2 <- ergm(
                ergm_model_2_params,
                # ERGMs with this many parameters can take a moment to compute
                control = control.ergm(
                  MCMC.samplesize = 5000, 
                  MCMC.interval = 8000, 
                  MCMC.burnin = 1 * 1000 * 1000, 
                  main.method = "Stochastic-Approximation",
                  
                  parallel = 8,           
                  parallel.type = "PSOCK"
                )
)
rough_estimates <- coef(noisy_model_2)
model_2 <- ergm(
  ergm_model_2_params,
  control = control.ergm(
    MCMC.samplesize = 5000, 
    MCMC.interval = 10 * 1000,
    MCMC.burnin = 100 * 1000, 
    init = rough_estimates,
    
    parallel = 8,           
    parallel.type = "PSOCK"
  )
)

# View the results to check which hypotheses are supported
summary(model_2)

saveRDS(model_2, file = "model_2_experimental.rds")

# model_2 <- readRDS(file = "model_2.rds")

library(coda)

# 2. Your model saves the data from all 8 cores in a list called "sample". 
# We use [[1]] to extract only the simulation data from the very first core.
single_thread_data <- model_2$sample[[1]]

# 3. Plot the trace and density for this single thread
plot(single_thread_data)

model_2_gof <-  gof(
  model_2,
  GOF = ~ idegree + odegree + espartners + distance
)
model_2_gof

ergm:::plot.gof(model_2_gof)
