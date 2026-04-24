library(ggplot2)
library(sf)
library(cowplot)
library(stringr)
library(dplyr)
library(tidyr)
library(xtable)
library(kableExtra)
map_probability <- readRDS("ProcessedData/map_probability.rds")



compute_indirect_risk <- function(
    trans_mat,
    county_name = "GAINES",
    threshold = .5
) {

  j <- toupper(county_name)  # source county j

  if (!(j %in% colnames(trans_mat))) stop("Source county not found in colnames(trans_mat).")

  i_names <- rownames(trans_mat)
  
  risk_table <- data.frame(
    County = i_names,
    first_gen_risk = NA_real_,
    second_gen_risk = NA_real_,
    stringsAsFactors = FALSE
  )
  
  P_ij <- trans_mat[ ,j]   # j -> i

  k_all <- names(P_ij)[!is.na(P_ij) & P_ij >= threshold & names(P_ij) != j]

  P_sg <- rep(NA_real_, length(P_ij))
  names(P_sg) <- i_names



  for (i in i_names) {
    # enforce k != i from paper equation
    k_use <- setdiff(k_all, i)

    if (length(k_use) == 0) {
      p2 <- 0
    } else {
      vals <- trans_mat[i, k_use] * trans_mat[k_use, j]  # P_ik * P_kj
      p2   <- mean(vals, na.rm = TRUE)
    }
    
    p1 <- P_ij[i]
    P_sg[i] <- min(1, p1 + p2) # combined first + second gen
    
    risk_table[risk_table$County == i, "first_gen_risk"] <- p1
    risk_table[risk_table$County == i, "second_gen_risk"] <- P_sg[i]
    
  }

  P_sg[j] <- NA_real_  # self county not shown
  
  
  risk_table <- risk_table %>%
    arrange(County)
  
  print(xtable(risk_table), 
        include.rownames = FALSE,
        sanitize.text.function = identity)
  return(P_sg)

  }



############### 

figure01 <- function(
    method = 7, 
    strategy = 0, 
    county = "Gaines",
    threshold = 0.5, 
    map_data = map_probability,
    out_dir = "Figures/"
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0-0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[, county_upper]  #j -> i
  indirect_vec <- compute_indirect_risk(
    trans_mat = mat, 
    county_name = county_upper, 
    threshold = threshold
    )
  
  highlighted_geom <- map_data %>% 
    filter(
      toupper(County) == county_upper
      )
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(direct_vec[match(toupper(County), names(direct_vec))],
                       breaks = breaks, 
                       labels = bin_labels, 
                       include.lowest = TRUE,
                       right = FALSE),
      indirect_bin = cut(indirect_vec[match(toupper(County), names(indirect_vec))],
                         breaks = breaks,
                         labels = bin_labels, 
                         include.lowest = TRUE,
                         right = FALSE)
      )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(
      aes(fill = direct_bin),
      color = "gray40",
      size = 0.1
      ) +
    geom_sf(
      data = highlighted_geom, 
      fill = "blue",
      color = "black", 
      size = 0.3
      ) +
    scale_fill_manual(
      values = color_palette,
      drop = FALSE,
      name = "Outbreak Probability"
      ) +
    labs(
      title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.position = "none"  
    )
  
  # Plot B: Indirect
  p2 <- ggplot(map_data) +
    geom_sf(
      aes(fill = indirect_bin), 
      color = "gray40", 
      size = 0.1
      ) +
    geom_sf(
      data = highlighted_geom, 
      fill = "blue", 
      color = "black",
      size = 0.3
      ) +
    scale_fill_manual(
      values = color_palette,
      drop = FALSE,
      name = "Outbreak Probability"
      ) +
    labs(
      title = "B: Probability of First & Second-generation Outbreak"
      ) +
    theme_void() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  legend <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  row_plot <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir, "figure01_outbreak_1st_2ndGen.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure01()

######################




#######################

figure02 <- function(
    method = 7, 
    counties = c("Gaines"),
    strategies = c(1, 2, 3),
    map_data = map_probability,
    threshold = 0.5, 
    out_dir = "Figures/"
    ) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#1a9850")  # white to green
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  baseline_mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S0.rds"))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    baseline_vec <- compute_indirect_risk(baseline_mat, county_name = county_upper, threshold = threshold)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- compute_indirect_risk(mat, county_name = county_upper, threshold = threshold)
      reduction <- (baseline_vec - vec) * 100
      reduction[reduction < 0] <- 0
      
      map_data$reduction <- reduction[match(toupper(map_data$County), names(reduction))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = reduction), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white", high = "#1a9850", limits = c(0, max(reduction, na.rm = TRUE)),
          name = "Reduction in\nOutbreak Probability (%)", na.value = "gray80"
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    legend <- get_legend(
      maps[[1]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure02_reduction", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(maps), height = 6, dpi = 400)
  }
}


figure02()


######################


figure03 <- function(
    counties,
    method = 7,
    map_data = map_probability,
    strategy = 0, 
    threshold = 0.5,
    out_dir = "Figures/"
) {
  
  # Define breaks and labels
  breaks     <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  
  
  county_names <- toupper(map_data$County)
  
  mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", strategy, ".rds"))
  maps <- list()
  
  panel_labels <- LETTERS[1:length(counties)]
  
  for (i in seq_along(counties)) {
    county <- counties[i]
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    indirect_vec <- compute_indirect_risk(
      mat, 
      county_name = county_upper, 
      threshold = threshold
    )
    
    map_data$indirect_risk <- indirect_vec[match(county_names, names(indirect_vec))]
    
    # Cut values into bins
    map_data <- map_data %>%
      mutate(
        indirect_bin = cut(
          indirect_vec[match(toupper(County), names(indirect_vec))],
          breaks = breaks,
          labels = bin_labels,
          include.lowest = TRUE,
          right = FALSE
        )
      )
    
    p <- ggplot(map_data) +
      geom_sf(aes(fill = indirect_bin), color = "gray40", size = 0.1) +
      geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
      scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
      labs(
        title = paste0(panel_labels[i], ": ", str_to_title(county)),
        fill = "Outbreak Probability"
      ) +
      theme_minimal() +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none"
      )
    
    maps[[length(maps) + 1]] <- p
  }
  
  legend <- get_legend(
    maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
  )
  
  row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
  
  file_name <- paste0(out_dir, "figure_03_indirect_risk_multiple_counties", method,
                      "_S", strategy, "_", paste(tolower(counties), collapse = "_"), ".png")
  ggsave(file_name, final_plot, width = 6 * length(counties), height = 6, dpi = 400)
  
}



figure03(counties = c("POLK", "MONTAGUE", "LIMESTONE"))
figure03(c("KING", "HALL", "THROCKMORTON"))


################

figure_S01 <- function(
    map_data = map_probability,
    out_dir = "Figures/"
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_mmr <- c("#d73027", "#fee08b", "#1a9850")     # for baseline MMR
  color_increase <- c("white", "#1a9850")             # for MMR increases
  
  # Compute increases (clamped)
  map_data <- map_data %>%
    mutate(
      delta_mmr1 = (MMR1 - MMR) * 100,
      delta_mmr2 = (MMR2 - MMR) * 100,
      delta_mmr3 = (MMR3 - MMR) * 100
    )
  
  plots <- list(
    ggplot(map_data) +
      geom_sf(aes(fill = MMR), color = "gray40", size = 0.1) +
      scale_fill_gradientn(
        colors = color_mmr,
        limits = c(0.8, 1),
        name = "MMR"
      ) +
      labs(title = "A: Baseline MMR (%)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr1), color = "gray40", size = 0.1) +
      scale_fill_gradientn(
        colors = color_increase, 
        limits = c(0, 20), 
        name = "Increase in MMR (%)"
      ) +
      labs(title = "B: MMR Vaccine Rate Increase by Strategy 01") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr2), color = "gray40", size = 0.1) +
      scale_fill_gradientn(
        colors = color_increase, 
        limits = c(0, 20), 
        name = "Increase in MMR (%)"
      ) +
      labs(title = "C: MMR Vaccine Rate Increase by Strategy 02") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr3), color = "gray40", size = 0.1) +
      scale_fill_gradientn(
        colors = color_increase, 
        #limits = c(0, 0.2),
        name = "Increase in MMR (%)"
      ) +
      labs(title = "D: MMR Vaccine Rate Increase by Strategy 03") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13))
  )
  
  grid <- plot_grid(plotlist = plots, ncol = 2, align = "hv")
  
  file_name <- paste0(out_dir, "Figure_S1_mmr_and_increase.png")
  ggsave(file_name, grid, width = 12, height = 10, dpi = 400)
  
  return(grid)
}



figure_S01()


##################

figure_S02 <- function(
    method = 7, 
    counties = c("Gaines"), 
    strategies = c(1, 2, 3),
    map_data = map_probability, 
    threshold = 0.5, 
    out_dir = "Figures/"
    ) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Discrete bins
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  names(color_palette) <- bin_labels
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  county_names <- toupper(map_data$County)
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- mats[[i]]
      indirect_vec <- compute_indirect_risk(
        mat, 
        county_name = county_upper, 
        threshold = threshold
        )
      
      bin_vec <- cut(indirect_vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
      bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(indirect_vec))], levels = bin_labels)
      map_data$bin <- bin_fact
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = bin), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_manual(
          values = color_palette,
          name = "Outbreak Probability",
          drop = FALSE
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    legend <- get_legend(
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure_S02_outbreak_2ndGen_strategy", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
  }
}



figure_S02()

#################



