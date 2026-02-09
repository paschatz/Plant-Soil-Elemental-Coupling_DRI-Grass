# Clean environment:
rm(list=ls())
getOption("device")
sessionInfo()

# Libraries:
library(tidyverse)
library(gridExtra)
library(nlme)
library(viridis)
library(cowplot)
library(grid)
library(Hmisc)
library(EcoCoupleR)
library(patchwork)
library(ggrepel)
library(georefdatar)
library(emmeans)

# Load data:
dat <- read_csv("data/DRIGrass_tidy.data.csv", show_col_types = FALSE)

# Functions:
cvar <- function(x) (sd(x)/mean(x))

se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# Store soil and plant variables:
soils.var <- names(dat) %>%
  keep(~ str_detect(.x, "Soil"))

plant.var <- names(dat) %>%
  keep(~ str_detect(.x, "Plant"))

# Set micro and macronutrients:
soil.macro <- c("Soil_N", "Soil_P", "Soil_K", "Soil_Ca", "Soil_Mg", "Soil_S")

soil.micro <- c("Soil_Fe", "Soil_B", "Soil_Mn", "Soil_Zn", "Soil_Cu")

plant.macro <- c("Plant_C", "Plant_N", "Plant_P", "Plant_K", "Plant_Ca", "Plant_Mg", "Plant_S", "Plant_Si")

plant.micro <- c("Plant_Cl", "Plant_Fe", "Plant_Mn", "Plant_Na", "Plant_Zn", "Plant_Cu")

pl.macro.names <- c("Carbon", "Nitrogen", "Phosphorus", "Potassium", "Calcium", "Magnesium", "Sulphur", "Silicon")

pl.micro.names <- c("Chlorine", "Iron", "Manganese", "Sodium", "Zinc" , "Copper")

# -------------------------------------------------------------------------------------- #
# Proceed with coupling
# -------------------------------------------------------------------------------------- #
data_coupling <- dat

cor.method <- "spearman"

# Function to calculate coupling per treatments (soil and plants):
compute_coupling_summary <- function(data, iteration) {
  map_dfr(unique(data$WHS), function(whs) {
    
    dat_soil <- data %>%
      filter(WHS == whs) %>%
      select(all_of(soils.var)) %>%
      as.matrix()
    
    dat_plant <- data %>%
      filter(WHS == whs) %>%
      select(all_of(plant.var)) %>%
      as.matrix()
    
    cor_soil <- rcorr(dat_soil, type = cor.method)
    cor_plant <- rcorr(dat_plant, type = cor.method)
    
    coupling_soil <- eco_coupling(cor_soil$r, data_str = "matrix")$coupling
    coupling_plant <- eco_coupling(cor_plant$r, data_str = "matrix")$coupling
    
    tibble(
      WHS = whs,
      iteration = iteration,
      coupling_soil = coupling_soil,
      coupling_plant = coupling_plant
    )
  })
}

# Function to calculate coupling per element (soil and plants):
compute_coupling_elm <- function(data, iteration) {
  map_dfr(unique(data$WHS), function(whs) {
    
    dat_soil <- data %>%
      filter(WHS == whs) %>%
      select(all_of(soils.var)) %>%
      as.matrix()
    
    dat_plant <- data %>%
      filter(WHS == whs) %>%
      select(all_of(plant.var)) %>%
      as.matrix()
    
    cor_soil <- rcorr(dat_soil, type = cor.method)
    cor_plant <- rcorr(dat_plant, type = cor.method)
    
    diag(cor_soil$r) <- NA
    diag(cor_plant$r) <- NA
    
    soil_elm <- cor_soil$r %>%
      as.data.frame() %>%
      mutate(elm = rownames(.)) %>%
      summarise(across(-elm, ~ mean(abs(.), na.rm = TRUE))) %>%
      pivot_longer(cols = everything(), names_to = "elm", values_to = "elm_coupling") %>%
      mutate(compartment = "soil", WHS = whs, iteration = iteration)
    
    plant_elm <- cor_plant$r %>%
      as.data.frame() %>%
      mutate(elm = rownames(.)) %>%
      summarise(across(-elm, ~ mean(abs(.), na.rm = TRUE))) %>%
      pivot_longer(cols = everything(), names_to = "elm", values_to = "elm_coupling") %>%
      mutate(compartment = "plant", WHS = whs, iteration = iteration)
    
    bind_rows(soil_elm, plant_elm)
  })
}

# Observed values:
# Store the observed values which are itteration zero (0):
observed_summary <- compute_coupling_summary(data_coupling, 0)
observed_elm <- compute_coupling_elm(data_coupling, 0)

# Permutations
set.seed(123)
num_permutations <- 999

perm_trt <- map_dfr(1:num_permutations, function(i) {
  permuted_data <- data_coupling %>%
    group_by(WHS) %>%
    mutate(across(where(is.numeric), ~ sample(.x, replace = FALSE))) %>%
    ungroup()
  
  compute_coupling_summary(permuted_data, i)
})

perm_elm <- map_dfr(1:num_permutations, function(i) {
  permuted_data <- data_coupling %>%
    group_by(WHS) %>%
    mutate(across(where(is.numeric), ~ sample(.x, replace = FALSE))) %>%
    ungroup()
  
  compute_coupling_elm(permuted_data, i)
})

# Plot coupling per treatment:
coupling.trt <- perm_trt %>% 
  bind_rows(observed_summary) %>%
  separate(WHS, into = c("Water", "Herbivores", "Year"), sep = "\\.", convert = TRUE) %>%
  # if herbivores 0 then -RH, if 1 then +RH
  mutate(treatment = case_when(
    Herbivores == 0 ~ paste(Water, "-RH", sep = " "),
    Herbivores == 1 ~ paste(Water, "+RH", sep = " "),
    TRUE ~ Water)) %>%
  mutate(Year = as.factor(Year),
         Herbivores = as.factor(Herbivores),
         Water = as.factor(Water),
         treatment = as.factor(treatment)) %>%
  # reorder based on treatment:
  mutate(treatment = factor(treatment, levels = c(
    "Control -RH", "Reduced -RH", "Altered Frequency -RH",
    "Control +RH", "Reduced +RH", "Altered Frequency +RH")))

# Reorder treatments for custom legend layout
coupling.trt <- coupling.trt %>%
  mutate(treatment = factor(
    treatment,
    levels = c(
      "Control -RH", 
      "Reduced -RH", 
      "Altered Frequency -RH", 
      "Control +RH", 
      "Reduced +RH", 
      "Altered Frequency +RH"
    )
  ))

# Separate observed and null
observed <- coupling.trt %>% 
  filter(iteration == 0)

null <- coupling.trt %>% 
  filter(iteration != 0)

# Compute null summaries
summary_df <- null %>%
  group_by(treatment, Year) %>%
  summarise(
    null_mean_soil = mean(coupling_soil),
    ci_low_soil = quantile(coupling_soil, 0.025),
    ci_high_soil = quantile(coupling_soil, 0.975),
    null_mean_plant = mean(coupling_plant),
    ci_low_plant = quantile(coupling_plant, 0.025),
    ci_high_plant = quantile(coupling_plant, 0.975),
    .groups = "drop") %>%
  left_join(
    observed %>% select(treatment, Year, obs_soil = coupling_soil, obs_plant = coupling_plant),
    by = c("treatment", "Year")) %>%
  mutate(
    sig_soil = ifelse(obs_soil < ci_low_soil | obs_soil > ci_high_soil, "*", ""),
    sig_plant = ifelse(obs_plant < ci_low_plant | obs_plant > ci_high_plant, "*", ""))

# Prepare long format for background jitter
null_long <- null %>%
  select(treatment, Year, coupling_soil, coupling_plant) %>%
  pivot_longer(cols = starts_with("coupling_"), names_to = "type", values_to = "value") %>%
  mutate(type = recode(type, coupling_soil = "Soil", coupling_plant = "Plant"))

# Pivot observed longer by putting together soil and plant coupling and creating a new variable (plant or soil)
observed_long <- observed %>%
  select(treatment, Year, coupling_soil, coupling_plant) %>%
  pivot_longer(cols = starts_with("coupling_"), names_to = "type", values_to = "value") %>%
  mutate(type = recode(type, coupling_soil = "Soil", coupling_plant = "Plant"))

################################################
# Figure 1 – Treatment Coupling (Soil and Plant)
################################################
# Shared scales and theme
y_axis <- scale_y_continuous(
  limits = c(0.25, 0.75),
  breaks = seq(0, 1, 0.25),
  expand = expansion(mult = c(0, 0.05)))

fill_scale <- scale_fill_viridis_d(option = "D")
color_scale <- scale_color_viridis_d(option = "D")

my_theme <- theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    strip.text = element_text(size = 11, face = "bold"),
    legend.spacing.y = unit(0.2, "cm")
  )


# Plot A – Soil
plot_soil <- ggplot() +
  geom_jitter(data = null_long %>% filter(type == "Soil"),
              aes(x = treatment, y = value, color = treatment),
              alpha = 0.2, width = 0.25, size = 0.7) +
  geom_errorbar(data = summary_df,
                aes(x = treatment, ymin = ci_low_soil, ymax = ci_high_soil),
                width = 0.2, color = "grey40") +
  geom_point(data = summary_df,
             aes(x = treatment, y = obs_soil, fill = treatment),
             shape = 21, size = 3, color = "black") +
  geom_text(data = summary_df,
            aes(x = treatment, y = obs_soil + 0.02, label = sig_soil),
            size = 5, color = "black") +
  facet_wrap(~Year, nrow = 1) +
  fill_scale + color_scale + y_axis + my_theme +
  labs(y = NULL, x = NULL, tag = "(a)", title = "Soil") +
  guides(color = "none", fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")

# Plot B – Plant
plot_plant <- ggplot() +
  geom_jitter(data = null_long %>% filter(type == "Plant"),
              aes(x = treatment, y = value, color = treatment),
              alpha = 0.2, width = 0.25, size = 0.7) +
  geom_errorbar(data = summary_df,
                aes(x = treatment, ymin = ci_low_plant, ymax = ci_high_plant),
                width = 0.2, color = "grey40") +
  geom_point(data = summary_df,
             aes(x = treatment, y = obs_plant, fill = treatment),
             shape = 21, size = 3, color = "black") +
  geom_text(data = summary_df,
            aes(x = treatment, y = obs_plant + 0.02, label = sig_plant),
            size = 5, color = "black") +
  facet_wrap(~Year, nrow = 1) +
  fill_scale + color_scale + y_axis + my_theme +
  labs(y = NULL, x = NULL, tag = "(b)", title = "Plant") +
  guides(color = "none", fill = "none") +  # remove second legend
  theme(legend.position = "none")

# Combine plots using patchwork and collect shared legend
panel_coupling <- wrap_plots(list(plot_soil, plot_plant), ncol = 2, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Convert to grob
panel_coupling_grob <- patchwork::patchworkGrob(panel_coupling)

# Add shared Y-axis label and display
fig1 <- gridExtra::grid.arrange(
  textGrob(
    expression("Global Elemental Coupling (|Spearman’s r|)"),
    rot = 90, gp = gpar(fontsize = 14), just = "centre"
  ),
  panel_coupling_grob,
  ncol = 2,
  widths = unit.c(unit(1.5, "cm"), unit(1, "null"))
)

# Show on screen
grid::grid.newpage()
grid::grid.draw(fig1)

# Save to file
ggsave("exports/DRIGrass_Fig1.tiff", plot = fig1, width = 9, height = 5, units = "in", dpi = 1200)

grid::grid.newpage()
grid::grid.draw(fig1)
dev.off()

# Elemental coupling ----
# Now proceed with elemental coupling for soils and plants:
coupling.elm <- perm_elm %>%
  bind_rows(observed_elm) %>%
  separate(WHS, into = c("Water", "Herbivores", "Year"), sep = "\\.", convert = TRUE) %>%
  # if herbivores 0 then -RH, if 1 then +RH
  mutate(treatment = case_when(
    Herbivores == 0 ~ paste(Water, "-RH", sep = " "),
    Herbivores == 1 ~ paste(Water, "+RH", sep = " "),
    TRUE ~ Water)) %>%
  mutate(Year = as.factor(Year),
         Herbivores = as.factor(Herbivores),
         Water = as.factor(Water),
         treatment = as.factor(treatment)) %>%
  # reorder based on treatment:
  mutate(treatment = factor(treatment, levels = c(
    "Control -RH", "Reduced -RH", "Altered Frequency -RH",
    "Control +RH", "Reduced +RH", "Altered Frequency +RH"
  )))

# 1. Separate observed and null
obs_elm <- coupling.elm %>% 
  filter(iteration == 0)

null_elm <- coupling.elm %>% 
  filter(iteration != 0)

# 2. Compute null summaries and merge with observed
summary_elm <- null_elm %>%
  group_by(treatment, Year, elm, compartment) %>%
  summarise(
    null_mean = mean(elm_coupling),
    ci_low = quantile(elm_coupling, 0.025),
    ci_high = quantile(elm_coupling, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    obs_elm %>%
      select(treatment, Year, elm, compartment, obs = elm_coupling),
    by = c("treatment", "Year", "elm", "compartment")
  ) %>%
  mutate(sig = ifelse(obs < ci_low | obs > ci_high, "*", ""))

######################################
# Figure 2 – Soil elemental coupling
######################################
# Use the raw variable names as default titles:
soil.vars.2 <- c("Soil_N", "Soil_P", "Soil_K", "Soil_Ca", "Soil_Mg", "Soil_S",   
                 "Soil_Fe", "Soil_B", "Soil_Mn", "Soil_Zn", "Soil_Cu")

soil_titles <- c("Nitrogen", "Phosphorus", "Potassium", "Calcium", "Magnesium", "Sulphur",
                 "Iron", "Boron", "Manganese", "Zinc", "Copper")

soil_tags <- paste0("(", letters[1:length(soil_titles)], ")")

# Create a look-up table for titles and tags
soil_title_map <- tibble(
  elm = soil.vars.2,
  title = soil_titles,
  tag = soil_tags)

# Filter soil-level coupling data
soil_data <- summary_elm %>% 
  filter(compartment == "soil") %>%
  semi_join(soil_title_map, by = "elm")

# Plot soil coupling function:
coupling_plot_soil <- function(data, title_map, show_legend = FALSE) {
  plots <- pmap(
    list(title_map$elm, title_map$title, title_map$tag, seq_along(title_map$elm)),
    function(elm_name, title, tag, idx) {
      dat <- data %>% filter(elm == elm_name)
      show_legend_now <- (idx == nrow(title_map)) && show_legend
      
      ggplot() +
        geom_jitter(data = null_elm %>%
                      filter(compartment == "soil", elm == elm_name) %>%
                      semi_join(dat, by = c("treatment", "Year", "elm")),
                    aes(x = treatment, y = elm_coupling, color = treatment),
                    alpha = 0.2, width = 0.2, size = 0.4, show.legend = FALSE) +
        geom_errorbar(data = dat,
                      aes(x = treatment, ymin = ci_low, ymax = ci_high),
                      width = 0.2, color = "grey40") +
        geom_point(data = dat,
                   aes(x = treatment, y = null_mean),
                   shape = 1, size = 2, color = "black") +
        geom_point(data = dat,
                   aes(x = treatment, y = obs, fill = treatment),
                   shape = 21, size = 3, color = "black", show.legend = show_legend_now) +
        geom_text(data = dat,
                  aes(x = treatment, y = obs + 0.02, label = sig),
                  size = 5, color = "black") +
        facet_wrap(~Year, ncol = 3) +
        labs(y = "", x = "", title = title, tag = tag) +
        scale_fill_viridis_d(option = "D", end = 0.9) +
        scale_color_viridis_d(option = "D", end = 0.9) +
        theme_bw(base_size = 13) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1.8, "lines"),
          legend.background = element_rect(fill = alpha("white", 0.9), colour = "transparent"),
          plot.title = element_text(hjust = 0.5),
          panel.spacing.x = unit(0, "lines"),
          plot.margin = margin(0, 0, 0, 0.2, "cm")
        )
    }
  )
  
  # Extract and remove legend
  legend_plot <- plots[[length(plots)]]
  legend <- get_legend(legend_plot)
  plots_clean <- map(plots, ~ . + theme(legend.position = "none"))
  
  # Arrange grid and draw
  arranged <- arrangeGrob(
    grobs = plots_clean, ncol = 3,
    left = textGrob("Soil elemental coupling (|Spearman’s r|)",
                    rot = 90, vjust = 1, gp = gpar(fontsize = 15))
  )
  
  ggdraw() +
    draw_grob(arranged) +
    draw_grob(legend, x = 0.75, y = 0.05, width = 0.12, height = 0.12)
}

# Generate the final figure
fig_soil <- coupling_plot_soil(
  soil_data,
  title_map = soil_title_map,
  show_legend = TRUE)

fig_soil

# Export soil coupling
ggsave("exports/fig2.tiff", fig_soil, width = 12, height = 15, dpi = 1200)

################################################
# Figure 4 & 5 – Elemental coupling (Plant)
################################################
null_elm <- coupling.elm %>% 
  filter(iteration != 0) %>%
  mutate(treatment = factor(
    treatment,
    levels = c(
      "Control -RH", 
      "Control +RH", 
      "Reduced -RH", 
      "Reduced +RH", 
      "Altered Frequency -RH",
      "Altered Frequency +RH"
    )
  ))

# Define plant macro/micro element names (as they appear in your data)
tags_macro <- paste0("(", letters[1:length(pl.macro.names)], ")")
tags_micro <- paste0("(", letters[1:length(pl.micro.names)], ")")

macro_data <- summary_elm %>% 
  filter(compartment == "plant", elm %in% plant.macro) %>%
  mutate(treatment = factor(
    treatment,
    levels = c(
      "Control -RH", 
      "Control +RH", 
      "Reduced -RH", 
      "Reduced +RH", 
      "Altered Frequency -RH",
      "Altered Frequency +RH"
    )
  ))

micro_data <- summary_elm %>% 
  filter(compartment == "plant", elm %in% plant.micro) %>%
  mutate(treatment = factor(
    treatment,
    levels = c(
      "Control -RH", 
      "Control +RH", 
      "Reduced -RH", 
      "Reduced +RH", 
      "Altered Frequency -RH",
      "Altered Frequency +RH"
    )
  ))

# Define plotting function
coupling_plot_plant_v2 <- function(data, title_list, tag_list, y_label) {
  elements <- unique(data$elm)
  
  # Generate individual plots
  plots <- map2(elements, seq_along(elements), ~ {
    dat <- data %>% filter(elm == .x)
    tag <- tag_list[.y]
    title <- title_list[.y]
    
    ggplot() +
      geom_jitter(data = null_elm %>%
                    filter(compartment == "plant", elm == .x) %>%
                    semi_join(dat, by = c("treatment", "Year", "elm")),
                  aes(x = treatment, y = elm_coupling, color = treatment),
                  alpha = 0.2, width = 0.2, size = 0.4) +
      geom_errorbar(data = dat,
                    aes(x = treatment, ymin = ci_low, ymax = ci_high),
                    width = 0.2, color = "grey40") +
      geom_point(data = dat,
                 aes(x = treatment, y = null_mean),
                 shape = 1, size = 2, color = "black") +
      geom_point(data = dat,
                 aes(x = treatment, y = obs, fill = treatment),
                 shape = 21, size = 3, color = "black") +
      geom_text(data = dat,
                aes(x = treatment, y = obs + 0.02, label = sig),
                size = 5, color = "black") +
      ylim(0, 0.8) +
      facet_wrap(~ Year, ncol = 3) +
      labs(y = "", x = "", title = title, tag = tag) +
      scale_fill_viridis_d(option = "D", end = 0.9) +
      scale_color_viridis_d(option = "D", end = 0.9) +
      theme_bw(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.8, "lines"),
        legend.background = element_rect(fill = alpha("white", 0.9), colour = "transparent"),
        plot.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        plot.margin = margin(0, 0, 0, 0.2, "cm")
      )
  })
  
  # Combine with patchwork and collect shared legend
  panel <- wrap_plots(plots, ncol = 2, guides = "collect") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Convert to grob
  panel_grob <- patchwork::patchworkGrob(panel)
  
  # Add shared Y-axis label
  fig <- gridExtra::grid.arrange(
    textGrob(
      y_label,
      rot = 90, gp = gpar(fontsize = 15), just = "centre"
    ),
    panel_grob,
    ncol = 2,
    widths = unit.c(unit(1.5, "cm"), unit(1, "null"))
  )
  
  return(fig)
}

# Figure 4:
fig_macro <- coupling_plot_plant_v2(
  macro_data,
  title_list = pl.macro.names,
  tag_list = tags_macro,
  y_label = expression("Plant elemental coupling (|Spearman’s r|)"))

# Figure 5:
fig_micro <- coupling_plot_plant_v2(
  micro_data,
  title_list = pl.micro.names,
  tag_list = tags_micro,
  y_label = expression("Plant elemental coupling (|Spearman’s r|)"))

# Export plant macronutrients
ggsave("exports/fig4.tiff", fig_macro, width = 12, height = 15, dpi = 1200)

# Export plant micronutrients
ggsave("exports/fig5.tiff", fig_micro, width = 12, height = 15, dpi = 1200)

# Linear mixed effect models:
aov(elm_coupling ~ Water * Herbivores * as.factor(Year), data = obs_elm %>% filter(compartment == "soil")) %>%
  anova() -> aov.soil

aov.soil

aov(elm_coupling ~ Water * Herbivores * as.factor(Year), data = obs_elm %>% filter(compartment == "plant")) %>%
  anova() -> aov.plant

aov.plant

# Elemental properties (atomic mass):
# Initialize empty list to store results
all_links <- list()

# Loop over treatments
for (trt in unique(data_coupling$WHS)) {
  
  dat_soil <- data_coupling %>%
    filter(WHS == trt) %>%
    select(all_of(soils.var)) %>%
    as.matrix()
  
  dat_plant <- data_coupling %>%
    filter(WHS == trt) %>%
    select(all_of(plant.var)) %>%
    as.matrix()
  
  cor_soil <- rcorr(dat_soil, type = "spearman")
  cor_plant <- rcorr(dat_plant, type = "spearman")
  
  diag(cor_soil$r) <- NA
  diag(cor_plant$r) <- NA
  
  cor_df_soil <- cor_soil$r %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
    filter(!is.na(weight)) %>%
    mutate(compartment = "soil", WHS = trt)
  
  cor_df_plant <- cor_plant$r %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
    filter(!is.na(weight)) %>%
    mutate(compartment = "plant", WHS = trt)
  
  # Store in list using the treatment as name
  all_links[[trt]] <- bind_rows(cor_df_soil, cor_df_plant)
}

# Function to extract atomic mass:
safe_aw <- possibly(
  function(sym) {
    # Get the atomic weight:
    result <- aw(sym)
    if (length(result) == 1 && is.numeric(result)) result else NA_real_
  },
  otherwise = NA_real_
)

# Combine all into one data frame
links_all_df <- bind_rows(all_links) %>%
  separate(WHS, into = c("Water", "Herbivores", "Year"), sep = "\\.", convert = TRUE) %>%
  # if herbivores 0 then -RH, if 1 then +RH
  mutate(treatment = case_when(
    Herbivores == 0 ~ paste(Water, "-RH", sep = " "),
    Herbivores == 1 ~ paste(Water, "+RH", sep = " "),
    TRUE ~ Water)) %>%
  # reorder based on treatment:
  mutate(treatment = factor(treatment, levels = c(
                      "Control -RH", "Reduced -RH", "Altered Frequency -RH",
                      "Control +RH", "Reduced +RH", "Altered Frequency +RH")),
         Year = as.factor(Year),
         Herbivores = as.factor(Herbivores),
         Water = as.factor(Water),
         from = as.factor(from),
         to = as.factor(to)) %>%
  mutate(weight = abs(weight)) %>%
  mutate(
    element = case_when(
      compartment == "soil" ~ str_split_fixed(from, "_", 2)[, 2],
      compartment == "plant" ~ {
        parts <- str_split(from, "_")
        map_chr(parts, ~ if (length(.x) == 3) .x[2] else tail(.x, 1))
      },
      TRUE ~ NA_character_)) %>%
  mutate(atomic_weight = map_dbl(element, safe_aw)) %>%
  mutate(log_weight = log(atomic_weight)) %>%
  # Make a new variable depending on if its macro or micronutrient:
  mutate(
    elm_type = case_when(
      from %in% plant.macro ~ "macronutrient",
      from %in% plant.micro ~ "micronutrient",
      from %in% soil.macro ~ "macronutrient",
      from %in% soil.micro ~ "micronutrient",
      TRUE ~ "other"))

lme(weight ~ log_weight * Water * Herbivores * Year, random = ~ 1 | Year,
    data = links_all_df %>% filter(compartment == "soil"),
    na.action = na.omit) -> lme_soil

anova(lme_soil)

summary.elm <- links_all_df %>%
  group_by(from, element, compartment, Water, atomic_weight, log_weight, elm_type) %>%
  summarise(elm.cor = mean(weight),
            se = se(weight)) %>%
  rowwise()

################################################
# Figure 3 – Atomic mass - soil elm
################################################

# Set theme:
theme_panel <- theme(panel.grid = element_blank(),
                     legend.title = element_blank(),
                     axis.text = element_text(size = 6),
                     axis.title = element_text(size = 8))
  
# Plot 1: Soil coupling vs atomic mass (Ambient)
p1 <- ggplot(summary.elm %>% filter(compartment == "soil", Water == "Control")) +
  geom_smooth(aes(x = log_weight, y = elm.cor),
              method = "lm", se = TRUE, color = "black", fill = "grey70") +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  labs(x = "Atomic mass (Da) - log",
       y = "Soil single-element coupling (|Spearman’s r|) – Ambient") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  theme_classic() +
  theme_panel

# Plot 2: Reduced amount
p2 <- ggplot(summary.elm %>% filter(compartment == "soil", Water == "Reduced")) +
  geom_smooth(aes(x = log_weight, y = elm.cor),
              method = "lm", se = TRUE, color = "black", fill = "grey70") +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  labs(x = "Atomic mass (Da) - log",
       y = "Soil single-element coupling (|Spearman’s r|) – Reduced am.") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  theme_classic() +
  theme_panel

# Plot 3: Altered frequency
p3 <- ggplot(summary.elm %>% filter(compartment == "soil", Water == "Altered Frequency")) +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  labs(x = "Atomic mass (Da) - log",
       y = "Soil single-element coupling (|Spearman’s r|) – Reduced freq.") +
  theme_classic() +
  theme_panel

# Proceed with coupling comparisons between treatments:
# Plot 4: Control vs Reduced amount
df_wide.fig.1 <- summary.elm %>%
  filter(compartment == "soil", Water %in% c("Control", "Reduced")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p4 <- ggplot(df_wide.fig.1, aes(x = cor_Control, y = cor_Reduced, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Soil single-element coupling (|Spearman’s r|) – Ambient", y = "Soil single-element coupling (|Spearman’s r|) – Reduced am.") +
  xlim(0.30, 0.55) + 
  ylim(0.30, 0.55) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

# Plot 5: Control vs Altered frequency
df_wide.fig.2 <- summary.elm %>%
  filter(compartment == "soil", Water %in% c("Control", "Altered Frequency")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p5 <- ggplot(df_wide.fig.2, aes(x = cor_Control, y = `cor_Altered Frequency`, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Soil single-element coupling (|Spearman’s r|) – Ambient", y = "Soil single-element coupling (|Spearman’s r|) – Reduced freq.") +
  xlim(0.30, 0.55) + 
  ylim(0.30, 0.55) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

# Plot 6: Reduced amount vs Altered frequency
df_wide.fig.3 <- summary.elm %>%
  filter(compartment == "soil", Water %in% c("Reduced", "Altered Frequency")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p6 <- ggplot(df_wide.fig.3, aes(x = cor_Reduced, y = `cor_Altered Frequency`, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Soil single-element coupling (|Spearman’s r|) – Reduced am.", y = "Soil single-element coupling (|Spearman’s r|) – Reduced freq.") +
  xlim(0.30, 0.55) + 
  ylim(0.30, 0.55) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

final_panel <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(5, 5, 5, 5))

# View it
print(final_panel)

ggsave("exports/fig3.tiff", final_panel, width = 10.5, height = 7.5, dpi = 1200)

################################################
# Figure S5 – Atomic mass - plant elm
################################################

# Repeat for plants:

# Plot 1: Soil coupling vs atomic mass (Ambient)
p1.pl <- ggplot(summary.elm %>% filter(compartment == "plant", Water == "Control")) +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  labs(x = "Atomic mass (Da) - log",
       y = "Plant single-element coupling (|Spearman’s r|) – Ambient") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  theme_classic() +
  theme_panel

# Plot 2: Reduced amount
p2.pl <- ggplot(summary.elm %>% filter(compartment == "plant", Water == "Reduced")) +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  labs(x = "Atomic mass (Da) - log",
       y = "Plant single-element coupling (|Spearman’s r|) – Reduced am.") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  theme_classic() +
  theme_panel

# Plot 3: Altered frequency
p3.pl <- ggplot(summary.elm %>% filter(compartment == "plant", Water == "Altered Frequency")) +
  geom_errorbar(aes(x = log_weight, ymin = elm.cor - se, ymax = elm.cor + se),
                width = 0.08, color = "darkgrey") +
  geom_point(aes(x = log_weight, y = elm.cor), size = 2, color = "black") +
  geom_text_repel(aes(x = log_weight, y = elm.cor, label = element),
                  size = 3, fontface = "italic", max.overlaps = 100,
                  nudge_x = 0.05, nudge_y = 0.01,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey30") +
  xlim(2.3, 4.3) + 
  ylim(0.30, 0.60) +
  labs(x = "Atomic mass (Da) - log",
       y = "Plant single-element coupling (|Spearman’s r|) – Reduced freq.") +
  theme_classic() +
  theme_panel

# Proceed with coupling comparisons between treatments:
# Plot 4: Control vs Reduced amount
df_wide.fig.1.plant <- summary.elm %>%
  filter(compartment == "plant", Water %in% c("Control", "Reduced")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p4.pl <- ggplot(df_wide.fig.1.plant, aes(x = cor_Control, y = cor_Reduced, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Plant single-element coupling (|Spearman’s r|) – Ambient", y = "Plant single-element coupling (|Spearman’s r|) – Reduced am.") +
  xlim(0.30, 0.6) + 
  ylim(0.30, 0.6) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

# Plot 5: Control vs Altered frequency
df_wide.fig.2.plant <- summary.elm %>%
  filter(compartment == "plant", Water %in% c("Control", "Altered Frequency")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p5.pl <- ggplot(df_wide.fig.2.plant, aes(x = cor_Control, y = `cor_Altered Frequency`, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Plant single-element coupling (|Spearman’s r|) – Ambient", y = "Plant single-element coupling (|Spearman’s r|) – Reduced freq.") +
  xlim(0.30, 0.6) + 
  ylim(0.30, 0.6) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

# Plot 6: Reduced amount vs Altered frequency
df_wide.fig.3.plant <- summary.elm %>%
  filter(compartment == "plant", Water %in% c("Reduced", "Altered Frequency")) %>%
  select(element, atomic_weight, Water, elm.cor, elm_type) %>%
  pivot_wider(names_from = Water, values_from = elm.cor, names_prefix = "cor_")

p6.pl <- ggplot(df_wide.fig.3.plant, aes(x = cor_Reduced, y = `cor_Altered Frequency`, color = elm_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  geom_text_repel(aes(label = element), size = 3, fontface = "italic", max.overlaps = 100,
                  show.legend = FALSE) +
  labs(x = "Plant single-element coupling (|Spearman’s r|) – Reduced am.", y = "Plant single-element coupling (|Spearman’s r|) – Reduced freq.") +
  xlim(0.30, 0.6) + 
  ylim(0.30, 0.6) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme_panel

final_panel.plants <- (p1.pl | p2.pl | p3.pl) / (p4.pl | p5.pl | p6.pl) +
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(5, 5, 5, 5))

# View it
print(final_panel.plants)

ggsave("exports/fig6.tiff", final_panel.plants, width = 10.5, height = 7.5, dpi = 1200)

##############################################
# Linear mixed effect models for soil coupling:
##############################################
# Run mixed effect models for global coupling:
mixed.effect.df <- links_all_df %>%
     # make interaction term by putting together from-to
     mutate(interaction = paste(from, to, sep = "-"),
            second.int = paste(interaction, Water, sep = "-"),
            Year = as.factor(Year))

# Run overall model for all soil elements combined
total.soils <- lme(weight ~ Water * Herbivores * Year,
                   random = ~ 1 | interaction/second.int,
                   data = mixed.effect.df %>% filter(compartment == "soil"))

anova_total.soils <- anova(total.soils) %>%
  rownames_to_column("Effect") %>%
  as.data.frame()

# Now repeat per element:
# Prepare an empty list to store results
anova_results.soil <- list()
for (elm in unique(soils.var)) {
  df_sub <- mixed.effect.df %>%
    filter(compartment == "soil", from == elm)
  
  model <- lme(weight ~ Water * Herbivores * Year,
               random = ~ 1 | interaction/second.int,
               data = df_sub)
  
  anova_results.soil[[elm]] <- anova(model)
}

# Shared structure
effects <- rownames(anova_results.soil[[1]])
numDF <- anova_results.soil[[1]]$numDF
denDF <- anova_results.soil[[1]]$denDF

# Assemble per-element results
anova_df.soils <- tibble(Effect = effects, numDF = numDF, denDF = denDF)

for (elm in names(anova_results.soil)) {
  temp <- anova_results.soil[[elm]]
  anova_df.soils[[paste0(elm, "_F")]] <- temp$`F-value`
  anova_df.soils[[paste0(elm, "_p")]] <- temp$`p-value`
}

# Merge total and per-element tables
full_anova_df.soils <- left_join(anova_total.soils, anova_df.soils, by = "Effect")


# Export table:
write_csv(full_anova_df.soils, "exports/soil_element_anova_results.csv")

##############################################
# Linear mixed effect models for plant coupling:
##############################################
# Run mixed effect models for global coupling:

# Run overall model for all soil elements combined
total.plants <- lme(weight ~ Water * Herbivores * Year,
                   random = ~ 1 | interaction/second.int,
                   data = mixed.effect.df %>% filter(compartment == "plant"))

anova_total.plants <- anova(total.plants) %>%
  rownames_to_column("Effect") %>%
  as.data.frame()

# Now repeat per element:
# Prepare an empty list to store results
anova_results.plant <- list()

for (elm in unique(plant.var)) {
  
  df_sub <- mixed.effect.df %>%
    filter(compartment == "plant", from == elm)
  
  model <- lme(weight ~ Water * Herbivores * Year,
               random = ~ 1 | interaction/second.int,
               data = df_sub)
  
  anova_results.plant[[elm]] <- anova(model)
}

# Shared structure
effects <- rownames(anova_results.plant[[1]])
numDF <- anova_results.plant[[1]]$numDF
denDF <- anova_results.plant[[1]]$denDF

# Assemble per-element results
anova_df.plant <- tibble(Effect = effects, numDF = numDF, denDF = denDF)

for (elm in names(anova_results.plant)) {
  temp <- anova_results.plant[[elm]]
  anova_df.plant[[paste0(elm, "_F")]] <- temp$`F-value`
  anova_df.plant[[paste0(elm, "_p")]] <- temp$`p-value`
}

# Merge total and per-element tables
full_anova_df.plants <- left_join(anova_total.plants, anova_df.plant, by = "Effect")

# Export table:
write_csv(full_anova_df.plants, "exports/plant_element_anova_results.csv")

