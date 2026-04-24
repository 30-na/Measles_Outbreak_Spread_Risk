library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Load county-level map data
map_data <- readRDS("ProcessedData/map_county.rds")
texas_flows <- read_csv("ProcessedData/texas_county_flows_2019.csv", 
                        col_types = cols(.default = "c"))

# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
R0 <- 18


map_data <- map_data %>%
  mutate(
    MMR1 = ifelse(MMR < .90, .90, MMR),           # Strategy 1
    MMR2 = ifelse(MMR < .92, .92, MMR),           # Strategy 2
    MMR3 = pmin(MMR + .05, 1),                    # Strategy 3
    # MMR4 = ifelse(County == "GAINES", 0.90, MMR),         # Strategy 4
    # MMR5 = ifelse(County == "GAINES", 0.92, MMR),         # Strategy 5
    
    # Susceptible proportion P^S
    susceptible_prop = 1 - (efficacy * MMR),
    susceptible_prop1 = 1 - (efficacy * MMR1),
    susceptible_prop2 = 1 - (efficacy * MMR2),
    susceptible_prop3 = 1 - (efficacy * MMR3),
    
    
    # susceptible_pop_size = susceptible_prop * total,
    # susceptible_pop_size1 = susceptible_prop1 * total,
    # susceptible_pop_size2 = susceptible_prop2 * total,
    # susceptible_pop_size3 = susceptible_prop3 * total,
    
    # Major local outbreak probability (Anderson & Watson form) P^LO
    outbreak_prob = pmax(0, 1 - (1 / ((1 - efficacy * MMR) * R0))),
    outbreak_prob1 = pmax(0, 1 - (1 / ((1 - efficacy * MMR1) * R0))),
    outbreak_prob2 = pmax(0, 1 - (1 / ((1 - efficacy * MMR2) * R0))),
    outbreak_prob3 = pmax(0, 1 - (1 / ((1 - efficacy * MMR3) * R0)))
    
  )


# State-level average vaccination coverage for each strategy
state_avg_coverage <- map_data %>%
  summarise(
    avg_MMR = mean(MMR, na.rm = TRUE),
    avg_MMR1 = mean(MMR1, na.rm = TRUE),
    avg_MMR2 = mean(MMR2, na.rm = TRUE),
    avg_MMR3 = mean(MMR3, na.rm = TRUE)
  )

# Number of counties that would increase to at least 90% and 92%
counties_to_90 <- sum(map_data$MMR < 0.90, na.rm = TRUE)
counties_to_92 <- sum(map_data$MMR < 0.92, na.rm = TRUE)

cat("State-Level Average Vaccination Coverage:\n")
cat(sprintf("  Original MMR:  %.4f\n", state_avg_coverage$avg_MMR))
cat(sprintf("  Strategy 1 (≥90%%): %.4f\n", state_avg_coverage$avg_MMR1))
cat(sprintf("  Strategy 2 (≥92%%): %.4f\n", state_avg_coverage$avg_MMR2))
cat(sprintf("  Strategy 3 (+5%%):  %.4f\n\n", state_avg_coverage$avg_MMR3))

cat("Number of counties reaching threshold:\n")
cat(sprintf("  To ≥90%% (Strategy 1): %d counties\n", counties_to_90))
cat(sprintf("  To ≥92%% (Strategy 2): %d counties\n", counties_to_92))

# Internal Infection Proportion Functions (Final outbreak size P^I) ----


# PLOS
find_internal_infection_PLOS <- function(Vj, efficacy_rate = efficacy, R0_value = R0, tol = 1e-4) {
  if (is.na(Vj)) return(NA_real_)
  a <- 1 - efficacy_rate * Vj
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-a * R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}

# Interface (Lin. , Equation: 2.4)
find_internal_infection_lin <- function(v, R0_value = R0, tol = 1e-4, efficacy_rate = efficacy) {
  if (is.na(v)) return(NA_real_)
  a <- 1 - (efficacy_rate * v) 
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}


# Apply function to MMR values
map_data$internal_infection_prob <- sapply(map_data$MMR, find_internal_infection_lin)
map_data$internal_infection_prob1 <- sapply(map_data$MMR1, find_internal_infection_lin)
map_data$internal_infection_prob2 <- sapply(map_data$MMR2, find_internal_infection_lin)
map_data$internal_infection_prob3 <- sapply(map_data$MMR3, find_internal_infection_lin)



# Compute average outbreak size
# map_data <- map_data %>%
#   mutate(
#     avg_outbreak_size = internal_infection_prob * total,
#     avg_outbreak_size1 = internal_infection_prob1 * total,
#     avg_outbreak_size2 = internal_infection_prob2 * total,
#     avg_outbreak_size3 = internal_infection_prob3 * total
#   )

saveRDS(map_data, "ProcessedData/map_probability.rds")



############ Contact Matrix ###########
pop <- map_data$total  # county population
names(pop) <- map_data$County


contact_method7 <- function(flows_avg) {
  
  county_names <- sort(unique(c(flows_avg$county_o, flows_avg$county_d)))
  mat <- matrix(0,
                length(county_names),
                length(county_names),
                dimnames = list(county_names, county_names))
  
  for (k in seq_len(nrow(flows_avg))) {
    origin <- flows_avg$county_o[k]
    dest <- flows_avg$county_d[k]
    flow <- flows_avg$flow[k]
    
    if (origin %in% county_names && dest %in% county_names) {
      if (origin == dest) {
        mat[dest, origin] <- NA  # self-flows as NA
      } else {
        mat[dest, origin] <- flow
      }
    }
  }
  
  return(mat)
}

  
  

# Method7
flows_avg <- texas_flows %>%
  mutate(
    pop_flows = as.numeric(pop_flows),
    visitor_flows = as.numeric(visitor_flows),
    week_start = mdy(str_sub(date_range, 1, 8))
  ) %>%
  filter(
     week_start >= as.Date("2019-01-15"),
     week_start <= as.Date("2019-08-15")
  ) %>%
  group_by(
    geoid_o,
    geoid_d, 
    county_o, 
    county_d
  ) %>%
  summarize(
    flow = sum(visitor_flows, na.rm = TRUE),
    .groups = "drop"
  )

C7 <- contact_method7(flows_avg)


############ Transmission/ Outbreak Probability ############


compute_transmission_matrix <- function(Cij, map_data, strategy = 0, q = 0.9) {
  suffix <- ifelse(strategy == 0, "", as.character(strategy))
  
  county_names <- map_data$County
  psi <- map_data[[paste0("susceptible_prop", suffix)]]
  pmo <- map_data[[paste0("outbreak_prob", suffix)]]
  Pj <- map_data[[paste0("internal_infection_prob", suffix)]]
  
  n <- length(county_names)
  transmission_mat <- matrix(NA,
                             n,
                             n,
                             dimnames = list(county_names, county_names))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      name_i <- county_names[i]
      name_j <- county_names[j]
      
      base <- q * Pj[which(county_names == name_j)] *
        psi[which(county_names == name_i)] *
        pmo[which(county_names == name_i)]
      
        # j -> i
        # cij_val <- Cij[name_i, name_j]
      
        # i -> j
        cij_val <- Cij[name_j, name_i]
      
      if (is.na(cij_val)) {
        transmission_mat[name_i, name_j] <- NA
      } else {
        transmission_mat[name_i, name_j] <- 1 - (1 - base)^cij_val
      }
      
    }
  }
  
  return(transmission_mat)
}


# Method 7 (Mobility Flows)
pij_M7_S0 <- compute_transmission_matrix(C7, map_data)
pij_M7_S1 <- compute_transmission_matrix(C7, map_data, strategy = 1)
pij_M7_S2 <- compute_transmission_matrix(C7, map_data, strategy = 2)
pij_M7_S3 <- compute_transmission_matrix(C7, map_data, strategy = 3)


saveRDS(pij_M7_S0, "ProcessedData/pij_M7_S0.rds")
saveRDS(pij_M7_S1, "ProcessedData/pij_M7_S1.rds")
saveRDS(pij_M7_S2, "ProcessedData/pij_M7_S2.rds")
saveRDS(pij_M7_S3, "ProcessedData/pij_M7_S3.rds")


