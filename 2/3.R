setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load the packages and source file
library(RSiena)
library(parallel)
source("printSiena.R") # Make sure this is in your working directory!

# 1. Data Loading ---------------------------------------------------------

# Adjacency matrices (No headers)
f1 <- as.matrix(read.csv("Glasgow/f1.csv", header = FALSE))
f2 <- as.matrix(read.csv("Glasgow/f2.csv", header = FALSE))
f3 <- as.matrix(read.csv("Glasgow/f3.csv", header = FALSE))

# Demographics and covariates
demographics <- read.csv("Glasgow/demographic.csv", header = TRUE)

# FIX: logdistance is a pure matrix just like the networks above
logdistance <- as.matrix(read.csv("Glasgow/logdistance.csv", header = FALSE))

# Alcohol has headers and an ID column (this is fine!)
alcohol <- as.matrix(read.csv("Glasgow/alcohol.csv", header = TRUE))

# 2. RSiena Data Processing -----------------------------------------------

# Dependent variable: Friendship network (129 nodes, 3 waves)
friendship <- as_dependent_rsiena(array(c(f1, f2, f3), dim = c(129, 129, 3)))

# Constant monadic covariates: Gender and Age
gender <- as_covariate_rsiena(demographics$gender, type = 'monadic')
age <- as_covariate_rsiena(demographics$age, type = 'monadic')

# Constant dyadic covariate: Log distance between houses
# (In RSiena wrappers, dyadic/matrix predictors are usually defined as oneMode)
distance <- as_covariate_rsiena(logdistance, type = 'oneMode')

# 3. Handle Alcohol -------------------------------------------------------
# Varying monadic covariate: Alcohol (3 waves)
alcohol_covar <- as_covariate_rsiena(alcohol, type = 'monadic')

# 4. Combine into RSiena Data Object --------------------------------------
mydata <- make_data_rsiena(friendship, gender, age, distance, alcohol_covar)
mydata

# 5. Generate Initial Report (Jaccard Index) ------------------------------
calc_jaccard <- function(net1, net2) {
  # Find ties that exist in both (maintained)
  maintained <- sum(net1 == 1 & net2 == 1, na.rm = TRUE)
  # Find ties that exist in at least one (total potential)
  total_ties <- sum(net1 == 1 | net2 == 1, na.rm = TRUE)
  
  return(maintained / total_ties)
}

# Print the Jaccard indices directly to the console
cat("Jaccard Period 1 to 2:", calc_jaccard(f1, f2), "\n")
cat("Jaccard Period 2 to 3:", calc_jaccard(f2, f3), "\n")

# 1. Create the initial effects object based on your data
myeff <- getEffects(mydata)

# 2. Add Endogenous Controls (Network Structure)
# Note: outdegree (density) and recip (reciprocity) are already included by default
myeff <- includeEffects(myeff, transTrip, cycle3)

# 3. Add Exogenous Controls (Demographics)
# Controlling for gender dynamics
myeff <- includeEffects(myeff, egoX, altX, sameX, interaction1 = "gender")

# Controlling for base alcohol tendencies (sender and receiver effects)
myeff <- includeEffects(myeff, egoX, altX, interaction1 = "alcohol_covar")
# (Optional) You can also add age controls if desired:
# myeff <- includeEffects(myeff, simX, interaction1 = "age")

# 4. Add Hypothesis Tests!
# H(i): Popularity
myeff <- includeEffects(myeff, inPop)

# H(ii): Alcohol Similarity
myeff <- includeEffects(myeff, simX, interaction1 = "alcohol_covar")

# H(iii): Distance / Proximity
myeff <- includeEffects(myeff, X, interaction1 = "distance")

# 5. Review your selected model specification
myeff

## Model estimation --------------------------------------------------------
# 5. Specify the parameters of the algorithm
myAlgorithm <- set_algorithm_saom(seed=1908, n3=3000, nsub=4)

# 6. Estimate the model
model0 <- siena07(myAlgorithm,
                  data = mydata,
                  effects = myeff,
                  returnDeps = TRUE,
                  batch = TRUE,
                  useCluster = TRUE,
                  nbrNodes = 8)

saveRDS(model0, file = "task3_model0.rds")
# printSiena(model0)

model0 <- readRDS("task3_model0.rds")
printSiena(model0)

# use parallel computation
cl <- makeCluster(4) # with 4 workers

# Indegree distribution
gofEvId <- test_gof(
  model0,
  verbose = FALSE,
  varName = "friendship", IndegreeDistribution,
  cluster = cl
)

# Outdegree distribution
gofEvOd <- test_gof(
  model0,
  verbose = FALSE,
  varName = "friendship", OutdegreeDistribution,
  cluster = cl
)

# Triad census
gofEvTC <- test_gof(
  model0,
  verbose = FALSE,
  varName = "friendship", TriadCensus,
  cluster = cl
)

# Geodesic distance
GeodesicDistribution <- function(
    i, data, sims, period, groupName,
    varName, levls = c(1:5, Inf), cumulative = TRUE
) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative) {
    gdi <- sapply(levls, function(i) {
      sum(a <= i)
    })
  }
  else {
    gdi <- sapply(levls, function(i) {
      sum(a == i)
    })
  }
  names(gdi) <- as.character(levls)
  gdi
}

gofEvGD <- test_gof(
  model0,
  verbose = FALSE,
  varName = "friendship", GeodesicDistribution
)

saveRDS(gofEvId, file = "task3_model0_gofEvId.rds")
saveRDS(gofEvOd, file = "task3_model0_gofEvOd.rds")
saveRDS(gofEvGD, file = "task3_model0_gofEvGD.rds")
saveRDS(gofEvTC, file = "task3_model0_gofEvTC.rds")

plot(gofEvId)
plot(gofEvOd)
plot(gofEvGD)
plot(gofEvTC, center = TRUE, scale = TRUE)

# 1. Redefine Alcohol as a dependent behavior variable
alcohol_beh <- as_dependent_rsiena(alcohol, type = "behavior")

# 2. Rebuild the RSiena Data Object
# Notice we drop the old 'alcohol_covar' and use 'alcohol_beh'
mydata <- make_data_rsiena(friendship, alcohol_beh, gender, age, distance)

myeff <- make_specification(mydata)

# =====================================================================
# SELECTION PART (Network Evolution)
# =====================================================================
# Endogenous Controls
myeff <- set_effect(myeff, list(transTrip, cycle3), depvar="friendship")

# Exogenous Controls (Gender)
myeff <- set_effect(myeff, list(egoX, altX, sameX), covar1="gender", depvar="friendship")

# Hypothesis (i): Popularity -> Friendship
myeff <- set_effect(myeff, list(inPop), depvar="friendship")

# Hypothesis (ii): Alcohol Similarity -> Friendship (Homophily)
# Note: covar1 is now targeting our new behavior variable "alcohol_beh"
myeff <- set_effect(myeff, list(simX, egoX, altX), covar1="alcohol_beh", depvar="friendship")

# Hypothesis (iii): Distance -> Friendship
myeff <- set_effect(myeff, list(X), covar1="distance", depvar="friendship")


# =====================================================================
# INFLUENCE PART (Behavior Evolution)
# =====================================================================
# Note: RSiena automatically includes the 'linear' and 'quad' (shape) 
# effects for behavior, which act as your basic endogenous controls.

# Exogenous Control: Gender's effect on alcohol consumption
# (It's good practice to control for basic demographics in behavior too)
myeff <- set_effect(myeff, list(effFrom), covar1="gender", depvar="alcohol_beh")

# Hypothesis (iv): Popularity -> Alcohol
# Does indegree in the friendship network predict alcohol behavior?
myeff <- set_effect(myeff, list(indeg), covar1="friendship", depvar="alcohol_beh")

# Hypothesis (v): Assimilation / Social Influence
# Do students adjust alcohol consumption to be similar to their friends?
myeff <- set_effect(myeff, list(avSim), covar1="friendship", depvar="alcohol_beh")

# Review the specification
myeff

# 1. Specify the algorithm (you can use your professor's wrapper here)
myAlgorithm <- set_algorithm_saom(seed=2024, n3=3000, nsub=4)

# 2. Estimate the model using 8 cores
modelCoev <- siena07(myAlgorithm, 
                     data = mydata, 
                     effects = myeff, 
                     returnDeps = TRUE,
                     batch = TRUE, 
                     useCluster = TRUE, 
                     nbrNodes = 8)
saveRDS(modelCoev, file = "task3_modelCoev.rds")

gofCoevId <- test_gof(
  modelCoev,
  verbose = FALSE,
  varName = "friendship", IndegreeDistribution,
  cluster = cl
)
plot(gofCoevId)
saveRDS(gofEvId, file = "task3_modelCoev_gofCoevId.rds")

# Outdegree distribution
gofCoevOd <- test_gof(
  modelCoev,
  verbose = FALSE,
  varName = "friendship", OutdegreeDistribution,
  cluster = cl
)
plot(gofCoevOd)
saveRDS(gofCoevOd, file = "task3_modelCoev_gofCoevOd.rds")

# Triad census
gofCoevTC <- test_gof(
  modelCoev,
  verbose = FALSE,
  varName = "friendship", TriadCensus,
  cluster = cl
)
plot(gofCoevTC)
saveRDS(gofCoevTC, file = "task3_modelCoev_gofCoevTC.rds")

# Geodesic distance
gofCoevGD <- test_gof(
  modelCoev,
  verbose = FALSE,
  varName = "friendship", GeodesicDistribution
)
plot(gofCoevGD)
saveRDS(gofCoevGD, file = "task3_modelCoev_gofCoevGD.rds")

# Behaviour distribution
gofCoevBeh <- test_gof(
  modelCoev,
  verbose = FALSE,
  varName = "alcohol_beh", BehaviorDistribution,
  cluster = cl
)
plot(gofCoevBeh)
saveRDS(gofCoevBeh, file = "task3_modelCoev_gofCoevBeh.rds")

stopCluster(cl)


