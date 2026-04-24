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

model1 <- ergm(glasgow_net1 ~ edges)
summary(model_1)

saveRDS(model_1, file = "task1_model1.rds")

# 4. Task ii: Fix decay to 0.3 and vary gwodegree coefficient over {-2, 2}

# --- Case 1: gwodegree coefficient = -2 ---
# We use offset() to fix the coefficient, and offset.coef to set its value.
# fixed = TRUE inside gwodegree() tells ergm that 0.3 is a fixed hyperparameter, 
# not a starting value for estimation.
model2_neg2 <- ergm(glasgow_net1 ~ edges + offset(gwodegree(0.3, fixed = TRUE)), 
                      offset.coef = -2)

print("Summary for gwodegree coef = -2")
summary(model2_neg2)

# Save the model
saveRDS(model2_neg2, file = "task1_model2_neg2.rds")


# --- Case 2: gwodegree coefficient = 2 ---
model2_pos2 <- ergm(glasgow_net1 ~ edges + offset(gwodegree(0.3, fixed = TRUE)), 
                      offset.coef = 2)

print("Summary for gwodegree coef = 2")
summary(model2_pos2)

# Save the model
saveRDS(model2_pos2, file = "task1_model2_pos2.rds")