set.seed(161)

library(sna)
library(network)

# Set working directory to the location of the current file (if using RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------------------------------
# 1. Define the Function
# ------------------------------------------------------------------------------------

global_nam_coefs <- data.frame()
global_nam_pvals <- data.frame()
global_nam_ses <- data.frame()

run_nam_for_classroom <- function(class_id, attr_raw, base_path = "./Lintner/") {
  file_w2 <- paste0(base_path, class_id, "_W2.csv")
  net_w2 <- as.matrix(read.csv2(file_w2, row.names = 1, check.names = FALSE))
  colnames(net_w2) <- gsub("^0", "", colnames(net_w2)) # Clean column names early
  
  w2_good_ids <- rownames(net_w2)[rowSums(is.na(net_w2)) < ncol(net_w2)]
  
  attr_base_ids <- attr_raw$studentID[
    attr_raw$classroomID == as.numeric(class_id) & 
      !is.na(attr_raw$literacy_end)
  ]
  ids_base <- intersect(attr_base_ids, w2_good_ids)
  
  attr_full_ids <- attr_raw$studentID[
    attr_raw$classroomID == as.numeric(class_id) & 
      !is.na(attr_raw$gender) & 
      !is.na(attr_raw$HISEI) & 
      !is.na(attr_raw$literacy_end)
  ]
  ids_full <- intersect(attr_full_ids, w2_good_ids)
  
  cat("\n======================================================\n")
  cat("             RESULTS FOR CLASSROOM:", class_id, "\n")
  cat("             Valid students (Base Model):", length(ids_base), "\n")
  cat("             Valid students (Full Model):", length(ids_full), "\n")
  cat("======================================================\n\n")
  
  # --- Setup & Run BASE Model ---
  attr_base <- attr_raw[attr_raw$studentID %in% ids_base, ]
  attr_base <- attr_base[match(ids_base, attr_base$studentID), ]
  attr_base$Intercept <- 1
  
  friendship_w2_base <- net_w2[as.character(ids_base), as.character(ids_base)]
  W_w2_base <- sweep(friendship_w2_base, 1, 
                     ifelse(rowSums(friendship_w2_base) == 0, 1, rowSums(friendship_w2_base)), 
                     FUN = "/")
  
  cat("--- Base NAM Model (Network Autocorrelation Only) ---\n")
  mod1 <- lnam(y = attr_base$literacy_end, x = as.matrix(attr_base$Intercept), W1 = W_w2_base)
  print(summary(mod1))
  
  attr_full <- attr_raw[attr_raw$studentID %in% ids_full, ]
  attr_full <- attr_full[match(ids_full, attr_full$studentID), ]
  attr_full$gender_num <- ifelse(attr_full$gender == "female", 1, 0)
  
  friendship_w2_full <- net_w2[as.character(ids_full), as.character(ids_full)]
  W_w2_full <- sweep(friendship_w2_full, 1, 
                     ifelse(rowSums(friendship_w2_full) == 0, 1, rowSums(friendship_w2_full)), 
                     FUN = "/")
  
  covars <- model.matrix(~ gender_num + HISEI, data = attr_full)
  
  cat("\n--- Full NAM Model (with Gender and HISEI) ---\n")
  mod2 <- lnam(y = attr_full$literacy_end, x = covars, W1 = W_w2_full)
  print(summary(mod2))
  cat("\n\n")
  
  # --- Extract Coefficients and P-values (From Full Model) ---
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
    Intercept = betas[1],
    Gender_Female = betas[2],
    HISEI = betas[3],
    Rho_Network = rho
  )
  
  row_pval <- data.frame(
    Classroom = class_id,
    Intercept = beta_pvals[1],
    Gender_Female_pval = beta_pvals[2],
    HISEI_pval = beta_pvals[3],
    Rho_Network_pval = rho_pval
  )
  
  row_se <- data.frame(
    Classroom = class_id,
    Intercept_se = beta_se[1],
    Gender_Female_se = beta_se[2],
    HISEI_se = beta_se[3],
    Rho_Network_se = rho_se
  )
  
  global_nam_coefs <<- rbind(global_nam_coefs, row_coef)
  global_nam_pvals <<- rbind(global_nam_pvals, row_pval)
  global_nam_ses <<- rbind(global_nam_ses, row_se) 
  
  return(list(
    classroom_id = class_id,
    n_students_base = length(ids_base),
    n_students_full = length(ids_full),
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
write.csv(global_nam_ses, "nam_results_se.csv", row.names = FALSE)

cat("Analysis complete! Results saved to nam_results_coefs.csv, nam_results_pvals.csv, and nam_results_se.csv.\n")