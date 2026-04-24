
# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
#(fuzzyjoin)
library(tigris)
library(sf)
library(stringdist)
library(purrr)
library(ggplot2)
library(readr)
#library(geosphere)


######### Replace mmr with average

# File paths
mmr_path_2024 <- "RawData/2023-2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx"
mmr_path_2023 <- "RawData/22-23-School-Vaccination-Coverage-by-District-and-County-K.xlsx"
mmr_path_2022 <- "RawData/2021-2022-School-Vaccination-Coverage-by-District-and-County-Kindergarten.xlsx"
mmr_path_2021 <- "RawData/2020-2021-School-Vaccination-Coverage-Levels-by-District-Private-School-and-County---Kindergarten.xlsx"
mmr_path_2020 <- "RawData/2019-2020-School-Vaccination-Coverage-Levels---Kindergarten.xlsx"

# Read sheets
sheet_2024 <- excel_sheets(mmr_path_2024)[2]
sheet_2023 <- excel_sheets(mmr_path_2023)[2]
sheet_2022 <- excel_sheets(mmr_path_2022)[2]
sheet_2021 <- excel_sheets(mmr_path_2021)[2]
sheet_2020 <- excel_sheets(mmr_path_2020)[2]

# Read each year and keep only County and MMR columns
mmr_2024 <- read_excel(mmr_path_2024, sheet = sheet_2024)
mmr_2023 <- read_excel(mmr_path_2023, sheet = sheet_2023, skip=2)
mmr_2022 <- read_excel(mmr_path_2022, sheet = sheet_2022, skip=2)
mmr_2021 <- read_excel(mmr_path_2021, sheet = sheet_2021, skip=2)
mmr_2020 <- read_excel(mmr_path_2020, sheet = sheet_2020, skip=2)

# Clean and select manually
mmr_2024 <- mmr_2024[, c("County", "MMR")] %>%
  mutate(
    MMR = as.numeric(MMR),
    County = toupper(County)
  )


mmr_2023 <- mmr_2023[, c("County", "MMR")] %>%
  mutate(
    MMR = as.numeric(MMR),
    County = toupper(County)
  )

mmr_2022 <- mmr_2022[, c("County", "MMR")] %>%
mutate(
  MMR = as.numeric(MMR),
  County = toupper(County)
)

mmr_2021 <- mmr_2021[, c("County", "MMR")] %>%
  mutate(
    MMR = as.numeric(MMR),
    County = toupper(County)
  )

mmr_2020 <- mmr_2020[, c("County", "MMR")] %>%
  mutate(
    MMR = as.numeric(MMR),
    County = toupper(County)
  )


# Rename MMR columns
names(mmr_2024)[2] <- "MMR_2024"
names(mmr_2023)[2] <- "MMR_2023"
names(mmr_2022)[2] <- "MMR_2022"
names(mmr_2021)[2] <- "MMR_2021"
names(mmr_2020)[2] <- "MMR_2020"

# Merge all years by County
df <- merge(mmr_2024, mmr_2023, by = "County", all = TRUE)
df <- merge(df, mmr_2022, by = "County", all = TRUE)
df <- merge(df, mmr_2021, by = "County", all = TRUE)
df <- merge(df, mmr_2020, by = "County", all = TRUE)



# Compute 5-year average (2020–2024)

df$MMR <- rowMeans(df[, c("MMR_2020", "MMR_2021","MMR_2022", "MMR_2023", "MMR_2024")], na.rm = TRUE)


mmr_avg <- df %>%
  select(County, MMR) %>%
  dplyr::filter(
    County != "TEXAS"
  )

# Save to CSV
write_csv(mmr_avg, "ProcessedData/county_mmr_5yr_average.csv")

mmr <- df %>%
  select(
    County,
    MMR
  )


tx_counties <- counties(state = "TX", year = 2024) %>%
  mutate(
    NAME = toupper(NAME)
  )


tx_counties_mmr <- left_join(
  tx_counties,
  mmr, 
  by = c("NAME" = "County")
) %>% 
  mutate(
    County = toupper(NAME)
  ) %>%
  select(
    County,
    MMR
  )



# Get the population
county_pop <- read_csv("RawData/counties_pop.txt") %>%
  mutate(
    FENAME = case_when(
      FENAME == "DE WITT" ~ "DEWITT",
      TRUE ~ FENAME)
  )


tx_counties_map <- tx_counties_mmr %>%
  left_join(
    county_pop,
    by = c("County" = "FENAME")
  ) %>%
  select(
    County,
    MMR,
    total
  )


saveRDS(tx_counties_map, "ProcessedData/map_county.rds")





###################### Population Flow
# Download code for the county2county weeky flow for 2019
# python download_weekly_data.py --start_year 2019 --start_month 1 --start_day 7 --end_year 2019 --end_month 12 --end_day 30 --output_folder weekly_flows  --county 



# Folder with weekly CSVs
file_list <- list.files("RawData/weekly_flows/county2county/", full.names = TRUE, pattern = "\\.csv$")

# Helper function to filter for Texas-only flows
is_texas_fips <- function(fips) {
  substr(fips, 1, 2) == "48"
}

# Read and filter all files
texas_flows <- lapply(file_list, function(file) {
  df <- read_csv(file, col_types = cols(.default = "c"))  # read as character to preserve leading 0s
  df <- df %>%
    filter(is_texas_fips(geoid_o) & is_texas_fips(geoid_d))
})

# Combine into one data frame
texas_flows_all <- bind_rows(texas_flows)

# Get county names from tigris
tx_counties <- counties(state = "TX", cb = TRUE, year = 2020) %>%
  st_drop_geometry() %>%
  select(GEOID, NAME)



# Add county names to origin and destination
texas_flows_named <- texas_flows_all %>%
  left_join(tx_counties, by = c("geoid_o" = "GEOID")) %>%
  rename(county_o = NAME) %>%
  left_join(tx_counties, by = c("geoid_d" = "GEOID")) %>%
  rename(county_d = NAME) %>%
  mutate(
    county_o = toupper(county_o),
    county_d = toupper(county_d)
  )

# Save final result with county names
write_csv(texas_flows_named, "ProcessedData/texas_county_flows_2019.csv")





