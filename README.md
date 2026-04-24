# Risk and Spatial Spread of a Measles Outbreak in Texas

This repository contains the code and data processing workflow for the manuscript **Risk and Spatial Spread of a Measles Outbreak in Texas**.

The project provides a county level modeling framework for estimating the risk and spatial spread of measles outbreaks across Texas. The model integrates county-level MMR vaccination coverage, county population data, and human mobility flows to estimate the probability that an outbreak in one county can seed major outbreaks in other counties. The framework is also used to evaluate vaccination improvement scenarios and alternative outbreak origin counties.

## Associated Manuscript

**Title:** Risk and Spatial Spread of a Measles Outbreak in Texas

**Authors:** Sina Mokhtar, Abhishek Pandey, Chad R. Wells, Martial L. Ndeffo-Mbah

**Repository:** [https://github.com/30-na/Measles_Outbreak_Spread_Risk](https://github.com/30-na/Measles_Outbreak_Spread_Risk)

## Overview

This repository supports a data driven analysis of measles outbreak risk in Texas. The main goals are to:

- Estimate the probability that a measles outbreak originating in one Texas county can seed major local outbreaks in other Texas counties.
- Evaluate first generation and second generation spatial spread of measles outbreaks across Texas counties.
- Compare model predictions with reported county level patterns from the 2025 West Texas measles outbreak.
- Assess the potential impact of MMR vaccination coverage under three vaccination scenarios.
- Simulate hypothetical outbreaks originating in selected low vaccination Texas counties.



## Data Sources

#### Vaccination Coverage
- **Source**: Texas Department of State Health Services  
- **Link**: [MMR Vaccination Coverage (DSHS)](https://www.dshs.texas.gov/immunizations/data/school/coverage)  
- File used: `2023–2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx`

#### County Population
- **Source**: Texas Legislative Council (2020 Census Redistricting Data)  
- **Link**: [Texas Capitol Data Portal](https://data.capitol.texas.gov/dataset/vtds)  
- File used: `Counties_Pop.txt`  
- Contains population data per county.

#### Mobility Flow Data
- **Source**: GeoDS COVID-19 US Flows Repository  
- **Link**: [https://github.com/GeoDS/COVID19USFlows](https://github.com/GeoDS/COVID19USFlows)  
- Weekly county-to-county flows were downloaded using:
 
```bash 

python download_weekly_data.py --start_year 2019 --start_month 1 --start_day 7 
--end_year 2019 \--end_month 12 --end_day 30 
--output_folder RawData/weekly_flows --county
```

#### Geometries Data
- **County Geometries**: TIGER/Line shapefiles via the `tigris` R package

## Model Description

To estimate the risk of a measles outbreak spreading in Texas, we considered the state as a network of counties linked through human mobility. The basis of the model is that the network of contacts between counties can provide estimates of the probability that a county outbreak can seed an outbreak in an interconnected county.

For each pair of counties, we calculated the probability that an outbreak could be seeded in county \( i \) given that an outbreak occurs in neighbouring county \( j \).

### Transmission Risk Model

To compute this probability, we first considered the probability of measles transmission through direct contact between an infectious and a susceptible individual as a set value \( q \).

The proportion of individuals in county \( j \) becomes infected during a local outbreak is denoted by \( P_j^I \), and the probability that an individual in county \( i \) is susceptible is denoted by \( P_i^S \). The probability that a single infected individual in county \( i \) causes a large outbreak in that county is denoted by \( P_i^{LO} \).

The probability that an outbreak in county \( j \) leads to an outbreak in county \( i \) through a single contact between the counties is given by:

```math
P_j^I P_i^S q P_i^{LO}
```

### Proportion Infected During a Local Outbreak

We assumed \( P_j^I \) to be equal to the proportion of individuals infected during the outbreak in county \( j \). This proportion is considered to be equal to the solution of the final outbreak size equation:

```math
P_j^I = (1 - M_j)(1 - e^{-R_0 P_j^I})
```

where \( R_0 \) is the measles basic reproductive number, and \( M_j \) is the immunity level in county \( j \).

In this analysis:

```math
R_0 = 18
```

### Susceptible Proportion

We assumed that \( P_i^S \) is equal to:

```math
P_i^S = 1 - M_i
```

We assumed that the measles immunity level of a given county \( M_i \) is a direct function of the proportion of vaccinated children younger than 10 years old.

We set:

```math
M_i = \epsilon V_i
```

where \( \epsilon \) is the MMR vaccine efficacy.

In this analysis, the MMR vaccine efficacy is:

```math
\epsilon = 0.97
```

### Probability of a Major Local Outbreak

We defined \( P_i^{LO} \) using the Anderson-Watson formula for the probability of a major outbreak of an SEIR-type pathogen in a susceptible population:

```math
P_i^{LO} = 1 - \frac{1}{(1 - M_i)R_0}
```


### First-Generation Outbreaks

If \( C_{ij} \) is the number of contact pairs that link county \( i \) and county \( j \), the probability that at least one contact pair causes a major outbreak in county \( i \) is given by:

```math
P_{ij} = 1 - (1 - P_j^I P_i^S q P_i^{LO})^{C_{ij}}
```

In this analysis:

```math
q = 0.9
```

We defined \( C_{ij} \) as the number of visits from county \( i \) to county \( j \) during 8 months, from mid-January to mid-August. To estimate \( C_{ij} \), we used between-county visitor flows in 2019 computed from anonymous mobile phone users.


### Second-Generation Outbreaks

A second-generation outbreak is a county outbreak resulting from contacts with infected individuals from a first-generation outbreak county.

If we denote by \( j \) the source county, the expected probability of a second-generation outbreak in county \( i \) is defined as:

```math
P_{ij}^{SG} = \frac{1}{n} \sum_{k \ne i} z_k P_{ik} P_{kj}
```

where \( k \) are potential first-generation outbreak counties.

```math
z_k = 1
```

if county \( k \) experienced a first-generation outbreak, and:

```math
z_k = 0
```

otherwise.

The expected probability value of a second-generation outbreak, \( P_{ij}^{SG} \), is defined as a conditional probability. We assumed that a first-generation outbreak occurred if:

```math
P_{kj} \geq 0.5
```

where \( n \) is the number of first-generation counties.

The analysis is limited to first-generation and second-generation outbreaks. Higher-order generations, such as third- or fourth-generation outbreaks, are not included.


## Vaccination Scenarios

The analysis evaluates one baseline scenario and three vaccination improvement strategies.

| Scenario | Description |
|---|---|
| Baseline | Uses observed county-level MMR vaccination coverage |
| Strategy 1 | Sets all counties with MMR coverage below 90% to exactly 90% |
| Strategy 2 | Sets all counties with MMR coverage below 92% to exactly 92% |
| Strategy 3 | Increases all county-level MMR coverage values by 5%, capped at 100% |


## Main Scripts

- `01_read_clean_merge_data.R`: Reads, cleans, and merges the vaccination, population, mobility, and county geometry data.

- `02_calculate_outbreak_probability.R`: Implements the transmission-risk model and calculates first-generation rutbreaks risk, and vaccination strategy.

- `03_plots.R`: Generates the maps, figures, and table used for the manuscript and supplementary materials.

## Run the Analysis

Run the scripts in order from the repository root directory.

```bash
Rscript scripts/01_read_clean_merge_data.R
Rscript scripts/02_calculate_outbreak_probability.R
Rscript scripts/03_plots.R
```


## Reproducibility Notes

- Raw data files should be placed in the appropriate `data/raw/` or `RawData/` directories before running the scripts.
- Large raw files not be included directly in this repository and need to be downloaded separately.
- Generated figures and processed outputs are saved under the `outputs/` directory.
- File paths may need to be adjusted depending on the local directory structure.

## Citation

If you use this repository, please cite the associated manuscript:

```bibtex
@article{mokhtar2026measles,
  title={Risk and Spatial Spread of a Measles Outbreak in Texas},
  author={Mokhtar, Sina and Pandey, Abhishek and Wells, Chad R. and Ndeffo-Mbah, Martial L.},
  year={2026},
  note={Manuscript in preparation}
}
```

This citation will be updated after publication.

## Contact

For questions about the model and manuscript, please contact the corresponding author:

**Martial L. Ndeffo-Mbah**  
Texas A&M University  
Email: m.ndeffo@tamu.edu

For questions about the repository or analysis, please contact:

**Sina Mokhtar**  
Department of Mathematics & Statistics  
Email: smokhtar@unm.edu



