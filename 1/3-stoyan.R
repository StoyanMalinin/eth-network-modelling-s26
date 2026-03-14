library(sna)
library(network)

# Set working directory to the location of the current file (if using RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------------------------------
# 1. Define the Function
# ------------------------------------------------------------------------------------

global_nam_coefs <- data.frame()
global_nam_pvals <- data.frame()

run_nam_for_classroom <- function(class_id, attr_raw, base_path = "./Lintner/") {
  file_w2 <- paste0(base_path, class_id, "_W2.csv")
  net_w2 <- as.matrix(read.csv2(file_w2, row.names = 1))
  
  attr_raw_good_ids <- attr_raw$studentID[
    attr_raw$classroomID == as.numeric(class_id) & 
      !is.na(attr_raw$gender) & 
      !is.na(attr_raw$HISEI) & 
      !is.na(attr_raw$literacy_end)
  ]
  
  w2_good_ids <- rownames(net_w2)[rowSums(is.na(net_w2)) < ncol(net_w2)]
  
  ids <- intersect(attr_raw_good_ids, w2_good_ids)
  attr <- attr_raw[attr_raw$studentID %in% ids, ]
  attr <- attr[match(ids, attr$studentID), ]
  
  if (class_id == "10" | class_id == "12") {
    colnames(net_w2) <- gsub("^X", "", colnames(net_w2))
  } else {
    colnames(net_w2) <- gsub("^X0", "", colnames(net_w2))  
  }
  
  friendship_w2 <- net_w2[as.character(ids), as.character(ids)]
  
  cat("\n======================================================\n")
  cat("             RESULTS FOR CLASSROOM:", class_id, "\n")
  cat("             Valid students analyzed:", length(ids), "\n")
  cat("======================================================\n\n")
  
  attr$gender_num <- ifelse(attr$gender == "female", 1, 0)
  
  W_w2 <- sweep(friendship_w2, 1, 
                ifelse(rowSums(friendship_w2) == 0, 1, rowSums(friendship_w2)), 
                FUN = "/")
  
  literacy_end <- attr$literacy_end
  
  cat("--- Base NAM Model (Network Autocorrelation Only) ---\n")
  mod1 <- lnam(y = literacy_end, W1 = W_w2)
  print(summary(mod1))
  
  covars <- model.matrix(~ gender_num + HISEI, data = attr)
  cat("\n--- Full NAM Model (with Gender and HISEI) ---\n")
  mod2 <- lnam(y = literacy_end, x = covars, W1 = W_w2)
  print(summary(mod2))
  cat("\n\n")
  
  # --- Extract Coefficients and P-values ---
  # The lnam object stores estimates and standard errors directly
  # We calculate the p-values manually using the Z-score
  
  # Exogenous variables (Intercept, gender_num, HISEI)
  betas <- mod2$beta
  beta_se <- mod2$beta.se
  beta_z <- betas / beta_se
  beta_pvals <- 2 * pnorm(-abs(beta_z))
  
  # Network autocorrelation (rho1)
  rho <- mod2$rho1
  rho_se <- mod2$rho1.se
  rho_z <- rho / rho_se
  rho_pval <- 2 * pnorm(-abs(rho_z))
  
  row_coef <- data.frame(
    Classroom = class_id,
    Gender_Female = betas[2],
    HISEI = betas[3],
    Rho_Network = rho
  )
  
  row_pval <- data.frame(
    Classroom = class_id,
    Gender_Female_pval = beta_pvals[2],
    HISEI_pval = beta_pvals[3],
    Rho_Network_pval = rho_pval
  )
  
  global_nam_coefs <<- rbind(global_nam_coefs, row_coef)
  global_nam_pvals <<- rbind(global_nam_pvals, row_pval)
  
  return(list(
    classroom_id = class_id,
    n_students = length(ids),
    model_base = mod1,
    model_full = mod2
  ))
}

# ------------------------------------------------------------------------------------
# 2. Main Script Execution
# ------------------------------------------------------------------------------------

attr_raw <- read.csv2("./Lintner/attr.csv", header = TRUE)
attr_raw[attr_raw == "" | attr_raw == " "] <- NA
attr_raw$literacy_end <- as.numeric(attr_raw$literacy_end)
attr_raw$HISEI <- as.numeric(attr_raw$HISEI)

all_classrooms <- list("01", "02", "03", "04", "05", "09", "10", "12")
all_nam_results <- lapply(all_classrooms, function(cid) {
  cat("Running NAM analysis for classroom:", cid, "...\n")
  run_nam_for_classroom(class_id = cid, attr_raw = attr_raw)
})

write.csv(global_nam_coefs, "nam_results_coefs.csv", row.names = FALSE)
write.csv(global_nam_pvals, "nam_results_pvals.csv", row.names = FALSE)

cat("Analysis complete! Results saved to nam_results_coefs.csv and nam_results_pvals.csv.\n")