# main_template.R

cat("Working directory:", getwd(), "\n")

output_dir <- file.path(getwd(), "outputs/script15")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Created directory:", output_dir, "\n")

source("DTvsfDT_alg.R")

rho <- 0.5
dt <- 15

keep_trav_cost <-c()

keep_rank <- list()
keep_seq <- list()

keep_margin <-c()



for(ii in 71:80){
  
  obs_gen_para <- read.csv(paste0("pattern/rho_", rho, "_", ii, ".csv"))
  the_result <- fDT_Alg_Save(obs_gen_para, dt)
  
  keep_trav_cost <- c(keep_trav_cost, the_result$trav_cost)
  
  keep_rank[[as.character(ii)]] <- the_result$rank_edges
  keep_seq[[as.character(ii)]] <- the_result$full_edge
  
  keep_margin <- c(keep_margin, the_result$avg_margin)
  
}

final_result <- list(trav_cost = keep_trav_cost,
                     rank = keep_rank, 
                     seq = keep_seq,
                     margin = keep_margin)

file_name <- file.path(output_dir, paste0("rho_", dt, "_", rho, ".rds"))
saveRDS(final_result, file = file_name)
