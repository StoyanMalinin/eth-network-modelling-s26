# Load the package
library(RSiena)

# 1. Data Loading ---------------------------------------------------------
# Be sure to set your working directory to the unzipped Glasgow folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the adjacency matrices (3 waves)
f1 <- as.matrix(read.csv("f1.csv", header = FALSE))
f2 <- as.matrix(read.csv("f2.csv", header = FALSE))
f3 <- as.matrix(read.csv("f3.csv", header = FALSE))

# Load demographic and covariate data
demographics <- read.csv("demographic.csv", header = TRUE)
logdistance <- as.matrix(read.csv("logdistance.csv", header = FALSE))

# Load alcohol data
alcohol <- as.matrix(read.csv("alcohol.csv", header = TRUE)) 

# 2. RSiena Data Processing -----------------------------------------------

# Dependent variable: Friendship network (129 nodes, 3 waves)
friendship <- sienaDependent(array(c(f1, f2, f3), dim = c(129, 129, 3)))

# Constant monadic covariates: Gender and Age
gender <- coCovar(demographics$gender)
age <- coCovar(demographics$age)

# Constant dyadic covariate: Log distance between houses
distance <- coDyadCovar(logdistance)

# Dependent behavior variable: Alcohol consumption (assumed for Co-evolution)
alcohol_beh <- sienaDependent(alcohol, type = "behavior")

# 3. Combine into RSiena Data Object --------------------------------------
mydata <- sienaDataCreate(friendship, alcohol_beh, gender, age, distance)
mydata

# 4. Generate Initial Report (includes Jaccard Index) ---------------------
# This will save a file named "Glasgow_Initial_Report.txt" to your folder
print01Report(mydata, modelname = "Glasgow_Initial_Report")