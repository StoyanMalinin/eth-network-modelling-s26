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
# 1. Define the Function
# ------------------------------------------------------------------------------------

global_coefs <- data.frame()
global_pvals <- data.frame()

run_qap_for_classroom <- function(class_id, attr_raw, base_path = "./Lintner/") {
  file_w1 <- paste0(base_path, class_id, "_W1.csv")
  file_w2 <- paste0(base_path, class_id, "_W2.csv")
  net_w1 <- as.matrix(read.csv2(file_w1, row.names = 1))
  net_w2 <- as.matrix(read.csv2(file_w2, row.names = 1))
  
  attr_raw_good_ids <- attr_raw$studentID[
    attr_raw$classroomID == as.numeric(class_id) & 
      !is.na(attr_raw$gender) & 
      !is.na(attr_raw$HISEI) & 
      !is.na(attr_raw$literacy_end)
  ]
  w1_good_ids <- rownames(net_w1)[rowSums(is.na(net_w1)) < ncol(net_w1)]
  w2_good_ids <- rownames(net_w2)[rowSums(is.na(net_w2)) < ncol(net_w2)]
  ids <- intersect(intersect(attr_raw_good_ids, w1_good_ids), w2_good_ids)
  
  attr <- attr_raw[attr_raw$studentID %in% ids, ]
  attr <- attr[match(ids, attr$studentID), ]
  
  if (class_id == "10" | class_id == "12") {
    colnames(net_w1) <- gsub("^X", "", colnames(net_w1))
    colnames(net_w2) <- gsub("^X", "", colnames(net_w2))
  } else {
    colnames(net_w1) <- gsub("^X0", "", colnames(net_w1))
    colnames(net_w2) <- gsub("^X0", "", colnames(net_w2))  
  }
  
  friendship_w1 <- net_w1[as.character(ids), as.character(ids)]
  friendship_w2 <- net_w2[as.character(ids), as.character(ids)]
  
  cat("\n======================================================\n")
  cat("             RESULTS FOR CLASSROOM:", class_id, "\n")
  cat("             Valid students analyzed:", length(ids), "\n")
  cat("======================================================\n\n")

  set.seed(161)
  permutations <- 5000
  reg1 <- netlogit(friendship_w2, friendship_w1, reps = permutations, nullhyp = "qapy")
  
  cat("--- Baseline Model (Wave 1) ---\n")
  print(summary(reg1))
  
  literacyReceiver <- matrix(attr$literacy_end, length(ids), length(ids), byrow = TRUE) # (i)
  sameGender <- outer(attr$gender, attr$gender, "==") * 1 # (ii)
  hiseiSender <- matrix(attr$HISEI, length(ids), length(ids), byrow = FALSE) # (iii)
  
  zm <- list(
    friendship_w1 = friendship_w1,
    literacyReceiver = literacyReceiver,
    sameGender = sameGender,
    hiseiSender = hiseiSender
  )
  
  reg2 <- netlogit(friendship_w2, zm, reps = permutations, nullhyp = "qapspp")
  
  cat("\n--- Full MR-QAP Model (Hypotheses i, ii, iii) ---\n")
  print(summary(reg2))
  cat("\n\n")
  
  # Store results for later
  
  ests <- reg2$coefficients
  # TODO (think about whether to use two-tailed or single-tailed p-value)
  pvals <- reg2$pgreqabs # This extracts the Pr(>=|b|) two-tailed p-value
  
  row_coef <- data.frame(
    Classroom = class_id,
    Friendship_from_W1 = ests[2],
    Literacy_Sender = ests[3],
    Gender_Match = ests[4],
    HISEI_Sender = ests[5]
  )
  row_pval <- data.frame(
    Classroom = class_id,
    Friendship_from_W1 = pvals[2],
    Literacy_Receiver = pvals[3],
    Gender_Match = pvals[4],
    HISEI_Sender = pvals[5]
  )
  
  global_coefs <<- rbind(global_coefs, row_coef)
  global_pvals <<- rbind(global_pvals, row_pval)
  
  return(list(
    classroom_id = class_id,
    n_students = length(ids),
    model_baseline = reg1,
    model_full = reg2
  ))
}

# ------------------------------------------------------------------------------------
# 2. Main Script Execution
# ------------------------------------------------------------------------------------

attr_raw <- read.csv2("./Lintner/attr.csv", header = TRUE)
attr_raw[attr_raw == "" | attr_raw == " "] <- NA  # Remember this step!
attr_raw$literacy_end <- as.numeric(attr_raw$literacy_end)

all_classrooms <- list("01", "02", "03", "04", "05", "09", "10", "12")
all_results <- lapply(all_classrooms, function(cid) {
  cat("Running analysis for classroom:", cid, "...\n")
  run_qap_for_classroom(class_id = cid, attr_raw = attr_raw)
})

write.csv(global_coefs, "qap_results_coefs.csv", row.names = FALSE)
write.csv(global_pvals, "qap_results_pvals.csv", row.names = FALSE)