# ------------------------------------------------------------------------------------
# Preparations
# ------------------------------------------------------------------------------------

# Load packages
# install.packages("sna")
# install.packages("network")
# install.packages("xtable")

library(sna)
library(network)
library(xtable)

# Set working directory to the location of the current file
# install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------------------------------
# Load data for classroom 10
# ------------------------------------------------------------------------------------

attr_raw <- read.csv2("./Lintner/attr.csv", header = TRUE)
attr_raw$literacy_end <- as.numeric(attr_raw$literacy_end)

net_w1 <- as.matrix(read.csv2("./Lintner/10_W1.csv", row.names = 1))
net_w2 <- as.matrix(read.csv2("./Lintner/10_W2.csv", row.names = 1))

attr_raw_good_ids <- attr_raw$studentID[attr_raw$classroomID == 10 & 
                                 !is.na(attr_raw$gender) & 
                                 !is.na(attr_raw$HISEI) & 
                                 !is.na(attr_raw$literacy_end)]

# IDs who provided network data (rows aren't all NA)
w1_good_ids <- rownames(net_w1)[rowSums(is.na(net_w1)) < ncol(net_w1)]
w2_good_ids <- rownames(net_w2)[rowSums(is.na(net_w2)) < ncol(net_w2)]

ids <- intersect(intersect(attr_raw_good_ids, w1_good_ids), w2_good_ids)

# Filter and Sort Attributes to match ids order
attr <- attr_raw[attr_raw$studentID %in% ids, ]
attr <- attr[match(ids, attr$studentID), ]

colnames(net_w1) <- gsub("^X", "", colnames(net_w1))
colnames(net_w2) <- gsub("^X", "", colnames(net_w2))
friendship_w1 <- net_w1[as.character(ids), as.character(ids)]
friendship_w2 <- net_w2[as.character(ids), as.character(ids)]

# ------------------------------------------------------------------------------------
# (1)
# ------------------------------------------------------------------------------------

set.seed(161)
permutations <- 5 * 1000
reg1 <- netlogit(friendship_w2, friendship_w1, reps = permutations, nullhyp = "qapy")

print(summary(reg1))

