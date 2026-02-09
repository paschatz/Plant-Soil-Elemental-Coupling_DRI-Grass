# Clean environemnt:
rm(list = ls())

# Load libraries:
library(tidyverse)
library(nlme)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(grid)
library(patchwork)
library(broom)

# Load dbroom# Load data:
datos <- read_csv("data/DRIGrass_tidy.data.csv")

# Store soil and plant variables:
soils.var <- names(datos) %>%
  keep(~ str_detect(.x, "Soil"))

plant.var <- names(datos) %>%
  keep(~ str_detect(.x, "Plant"))

# Set micro and macronutrients:
soil.macro <- c("Soil_N", "Soil_P", "Soil_K", "Soil_Ca", "Soil_Mg", "Soil_S")

soil.micro <- c("Soil_Fe", "Soil_B", "Soil_Mn", "Soil_Zn", "Soil_Cu")

plant.macro <- c("Plant_C", "Plant_N", "Plant_P", "Plant_K", "Plant_Ca", "Plant_Mg", "Plant_S", "Plant_Si")

plant.micro <- c("Plant_Cl", "Plant_Fe", "Plant_Mn", "Plant_Na", "Plant_Zn", "Plant_Cu")

pl.macro.names <- c("Carbon", "Nitrogen", "Phosphorus", "Potassium", "Calcium", "Magnesium", "Sulphur", "Silicon")

pl.micro.names <- c("Chlorine", "Iron", "Manganese", "Sodium", "Zinc" , "Copper")

# -------------------------------------------------------------------------------------- #
# Run multiple lme models for each variable (Sup. Table S1):
# -------------------------------------------------------------------------------------- #
# Initiate lists to store lme results:
multiple <- list()
thesecolumns <- c(soils.var, plant.var)

# Loop over columns to run lme models:
for(i in 1:length(thesecolumns)) {
  
  thiscolumn <- thesecolumns[i]
  datos$Y <- datos[[thiscolumn]]
  anovas <- lme(Y ~ Water.treatment * Year * Herbivores, 
                random = ~1|Plot.number, 
                data = datos, 
                na.action = na.omit)
  
  multiple[[i]] <- anovas
}

# Naming the list elements:
names(multiple) <- names(datos)[c(9:19,20:33)]

# Tidy equivalent of the base code
anova.table <- map2_dfc(multiple, names(multiple), ~ {
  out <- as.data.frame(anova(.x))
  out <- out[-1, ]  # remove intercept row
  out <- out[, c("numDF", "denDF", "F-value", "p-value")]
  colnames(out) <- paste0(.y, ".", colnames(out))  # e.g. Ca_PRS.F-value
  out}) %>%
  # round any numeric variable:
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  # Add rownames as a proper column
  rownames_to_column("Term") %>%
  # Replace the raw ANOVA rownames with clean names
  mutate(
    Term = recode(
      Term,
      "Water.treatment"                  = "Rainfall",
      "Year"                             = "Year",
      "Herbivores"                       = "Herbivores",
      "Water.treatment:Year"             = "Rainfall x Year",
      "Water.treatment:Herbivores"       = "Rainfall x Herbivores",
      "Year:Herbivores"                  = "Year x Herbivores",
      "Water.treatment:Year:Herbivores"  = "Rainfall x Year x Herbivores"
    )
  )

# Split into soil and plants:
soil.anova <- anova.table %>%
  select(Term, starts_with("Soil"))

plant.anova <- anova.table %>%
  select(Term, starts_with("Plant"))

# Export anova tables:
write.csv(soil.anova,"exports/DRIGrass_Table.S1.csv")

write.csv(plant.anova,"exports/DRIGrass_Table.S2.csv")
# -------------------------------------------------------------------------------------- #
##########################################################################################
# Figure S1: Soil nutrient bioavailability
##########################################################################################
# Desired order (all -RH first, then +RH)
wth_order <- c(
  "Ambient -RH", "Reduced amount -RH", "Reduced frequency -RH",
  "Ambient +RH", "Reduced amount +RH", "Reduced frequency +RH"
)

# Make WTH an ordered factor once (so axis *and* legend follow it)
plot.dat <- datos %>%
  mutate(WTH = factor(WTH, levels = wth_order))

# Plot soil elemental bioavailability:
element_plot <- function(var, title, tag, show_legend = FALSE) {
  ggplot(plot.dat, aes(x = WTH, y = .data[[var]], fill = WTH)) +
    geom_jitter(aes(colour = WTH), width = 0.15, alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, show.legend = show_legend) +
    geom_vline(xintercept = 3.5, color = "red", linetype = "dashed") +
    labs(y = "", x = "", title = title, tag = tag) +
    scale_fill_viridis_d(option = "D", end = 0.9) +  # fill = boxplot
    scale_colour_viridis_d(option = "D", end = 0.9) +  # colour = jitter
    facet_wrap(~Year, ncol = 5) +
    theme_bw(base_size = 13) +
    guides(
      fill   = guide_legend(nrow = 6, byrow = TRUE),
      colour = guide_legend(nrow = 6, byrow = TRUE)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 18, color = "black"),
      legend.key.size = unit(1.6, "lines"),
      legend.background = element_rect(fill = alpha("white", 0.9), colour = "transparent"),
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5),
      panel.spacing.x = unit(0, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.margin = margin(0, 0, 0, 0.2, "cm"),
      panel.border = element_rect(color = "black")
    )
}

# Use the raw variable names as default titles:
soil.vars.2 <- c("Soil_N", "Soil_P", "Soil_K", "Soil_Ca", "Soil_Mg", "Soil_S",   
                 "Soil_Fe", "Soil_B", "Soil_Mn", "Soil_Zn", "Soil_Cu")

soil_titles <- c("Nitrogen", "Phosphorus", "Potassium", "Calcium", "Magnesium", "Sulphur",
                 "Iron", "Boron", "Manganese", "Zinc", "Copper")

tags <- letters[1:length(soil.vars.2)]

# Create the list of plots: only the last one gets a legend
soil_plots <- map2(soil.vars.2, seq_along(soil.vars.2), ~ {
  show_legend <- (.y == length(soil.vars.2))  # show legend only in the last plot
  element_plot(.x, soil_titles[.y], paste0("(", tags[.y], ")"), show_legend = show_legend)
})

# Extract the legend from the last plot
legend_plot <- soil_plots[[length(soil_plots)]]
legend <- get_legend(legend_plot)

# Remove legends from all individual plots now (for safety)
soil_plots_clean <- map(soil_plots, ~ . + theme(legend.position = "none"))

# Arrange the plots
main_plot <- arrangeGrob(grobs = soil_plots_clean, ncol = 3,
                         left = textGrob(expression(paste("Soil nutrient biovailability (mg·m"^"-2", "·3 months"^"-1", ")")),
                                         rot = 90, vjust = 1, gp = gpar(fontsize = 15)))

# Combine with legend in bottom right corner
figs1 <- ggdraw() +
  draw_grob(main_plot) +
  draw_grob(legend, x = 0.78, y = 0.05, width = 0.15, height = 0.15)  # adjust x/y/size as needed

# Print
print(figs1)

# Export:
ggsave("exports/DRIGrass_Fig.S1.tiff", plot = figs1, width = 12, height = 15, dpi = 1200)

##########################################################################################
# Figure S2: Plant macronutrient bioavailability
##########################################################################################
# Plot function
element_plot_plants <- function(var, title, tag) {
  ggplot(plot.dat, aes(x = WTH, y = .data[[var]], fill = WTH)) +
    geom_jitter(aes(colour = WTH), width = 0.15, alpha = 0.5) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    geom_vline(xintercept = 3.5, color = "red", linetype = "dashed") +
    labs(y = NULL, x = NULL, title = title, tag = tag) +
    scale_fill_viridis_d(option = "D", end = 0.9, breaks = wth_order, drop = FALSE) +
    scale_colour_viridis_d(option = "D", end = 0.9, breaks = wth_order, drop = FALSE) +
    facet_wrap(~Year, ncol = 5) +
    theme_bw(base_size = 13) +
    guides(
      fill   = guide_legend(nrow = 2, byrow = TRUE),
      colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, color = "black"),
      legend.key.size = unit(1.6, "lines"),
      legend.box = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.9), colour = "transparent"),
      plot.title = element_text(hjust = 0.5),
      panel.spacing.x = unit(0, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      plot.margin = margin(0, 0, 0, 0.2, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
}

# Generate list of plots
tags_macro <- letters[1:length(plant.macro)]
plant_macro.plots <- map2(plant.macro, seq_along(plant.macro), ~ {
  element_plot_plants(.x, pl.macro.names[.y], paste0("(", tags_macro[.y], ")"))
})

# Combine with shared legend
panel_macro <- wrap_plots(plant_macro.plots, ncol = 2, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

# Convert patchwork object to grob
panel_macro_grob <- patchwork::patchworkGrob(panel_macro)

# Combine Y label and plot panel in one figure
figs2 <- gridExtra::grid.arrange(
  textGrob(
    expression("Plant elemental content (mg·kg"^"-1"*")"),
    rot = 90, gp = gpar(fontsize = 13), just = "centre"),
  panel_macro_grob,
  ncol = 2,
  widths = unit.c(unit(1.5, "cm"), unit(1, "null")))

grid.draw(figs2)

# Export
ggsave("exports/DRIGrass_Fig.S2.tiff",plot = figs2, width = 10, height = 15, units = "in", dpi = 1200)

dev.new()
grid.newpage()

##########################################################################################
# Figure S3: Plant micronutrient bioavailability
##########################################################################################
# Create micronutrient plots
tags_micro <- letters[1:length(plant.micro)]
plant_micro.plots <- map2(plant.micro, seq_along(plant.micro), ~ {
  element_plot_plants(.x, pl.micro.names[.y], paste0("(", tags_micro[.y], ")"))
})

# Create 2 empty placeholder plots
empty_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank())

# Add 2 blank plots to reach 8 total (4 rows x 2 columns)
plant_micro.filled <- c(plant_micro.plots, list(empty_plot, empty_plot))

# Combine plots with shared legend using patchwork
panel_micro <- wrap_plots(plant_micro.filled, ncol = 2, nrow = 4, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Convert patchwork object to grob
panel_micro_grob <- patchwork::patchworkGrob(panel_micro)

# Combine with shared Y-axis label
figs3 <- gridExtra::grid.arrange(
  textGrob(
    expression("Plant micronutrient bioavailability (mg·kg"^"-1"*")"),
    rot = 90, gp = gpar(fontsize = 16), just = "centre"
  ),
  panel_micro_grob,
  ncol = 2,
  widths = unit.c(unit(1.5, "cm"), unit(1, "null"))
)

grid.draw(figs3)

# Save to file
ggsave("exports/DRIGrass_Fig.S3.tiff", plot = figs3, width = 10, height = 15, units = "in", dpi = 1200)

dev.off()
grid.newpage()

################################################################################
# Proceed with Effect size:
# Define variables
# Set herbivores as factor
# Reference level is 0 (without herbivores)
datos$Herbivores <- factor(datos$Herbivores, levels = c(0, 1))

# Repeat for water
datos$Water.treatment <- factor(datos$Water.treatment, levels = c("Control", "Reduced", "Altered Frequency"))

# Store elements
elm.df <- c(soils.var, plant.var) # Ensure this list exists

# Initialize list to store results
results_list <- list()
k <- 1

# Loop oveer years:
for (year in unique(datos$Year)) {
  
  # Filter data for the specific year
  dat.year <- datos %>% 
    filter(Year == year)
    
  # Loop over elements:
  for (elm in elm.df) {
    f.water <- as.formula(paste0("scale(", elm, ") ~ Water.treatment"))
    f.herbivores <- as.formula(paste0("scale(", elm, ") ~ Herbivores"))
    
    # Run Model
    model.water <- lm(f.water, data = dat.year)
    mode.herbivores <- lm(f.herbivores, data = dat.year)
    
    # Extract Coefficients
    res.water <- tidy(model.water, conf.int = TRUE, conf.level = 0.95) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        Year = year,
        full_name = elm) %>%
      
      # Separate medium-element
      separate(
        full_name,
        into   = c("medium", "element"),
        sep    = "_",
        remove = TRUE)
    
    res.herbivores <- tidy(mode.herbivores, conf.int = TRUE, conf.level = 0.95) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        Year = year,
        full_name = elm) %>%
      # Separate medium-element
      separate(
        full_name,
        into   = c("medium", "element"),
        sep    = "_",
        remove = TRUE)
    
    res <- bind_rows(res.water, res.herbivores)
    
    results_list[[k]] <- res
    k <- k + 1
  }
}

# Combine results into one data frame
model_results <- bind_rows(results_list)

# Tidy up:
models.tidy <- model_results %>%
  mutate(p.value = round(p.value, 3),
         sig = ifelse(p.value < 0.05, "sig", "ns"))

# Prepare Data
plot_data <- models.tidy %>%
  mutate(
    # --- Label Logic ---
    Treatment_Label = case_when(
      term == "Water.treatmentReduced" ~ "Reduced Amount",
      term == "Water.treatmentAltered Frequency" ~ "Reduced Frequency",
      term == "Herbivores1" ~ "+RH", 
      TRUE ~ term
    ),
    
    # --- Panel Logic ---
    Panel = case_when(
      str_detect(term, "Water") ~ "Water Treatment",
      str_detect(term, "Herbivore") ~ "Herbivory",
      TRUE ~ "Other"
    ),
    
    # Force the order of elements
    element = factor(element, levels = sort(unique(element), decreasing = TRUE)),
    
    # Set significance
    sig_alpha = ifelse(p.value < 0.05, "Significant", "NS")) %>%
  # Filter 
  filter(Panel %in% c("Water Treatment", "Herbivory"))

# Plot
ggplot(plot_data, aes(x = estimate, y = element, color = Treatment_Label, group = Treatment_Label)) +
  
  # Reference Line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Error Bars & Points
  geom_pointrange(
    aes(xmin = conf.low, xmax = conf.high, alpha = sig_alpha),
    size = 0.2,
    position = position_dodge(width = 0.7),
    flatten = 0.1
  ) +
  
  # Faceting: 
  facet_grid(Panel + medium ~ Year, scales = "free_y") +
  
  # --- COLORS ---
  scale_color_manual(values = c(
    "Reduced Amount"    = "#D55E00",
    "Reduced Frequency" = "#0072B2",
    "+RH"= "#009E73")) +
  
  # Non-significant results are faded, Significant are solid
  scale_alpha_manual(values = c("NS" = 0.4, "Significant" = 1.0)) +
  
  theme_minimal() +
  labs(
    title = "",
    x = "Standardized Effect Size (ES)",
    y = NULL,
    color = "Treatment",
    alpha = "Significance"
  ) +
  
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    axis.title = element_text(size = 8, color = "black"),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text = element_text(color = "black", size = 6),
    legend.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 6, color = "black", face = "bold")
    
  )

# GGsave:
ggsave("exports/DRIGrass_Fig.SUP.tiff", width = 15, height = 15, units = "cm", dpi = 1200)

# Initialize list
results_list.global <- list()
k <- 1

# Loop over elements
for (elm in elm.df) {
  
  # DEFINING FORMULAS
  f.water <- as.formula(paste0("scale(", elm, ") ~ Water.treatment"))
  f.herbivores <- as.formula(paste0("scale(", elm, ") ~ Herbivores"))
  
  # RUN MODELS
  model.water <- lm(f.water, data = datos) 
  mode.herbivores <- lm(f.herbivores, data = datos)
  
  # TIDY WATER
  res.water <- tidy(model.water, conf.int = TRUE, conf.level = 0.95) %>%
    filter(term != "(Intercept)") %>%
    mutate(full_name = elm) %>%
    separate(full_name, into = c("medium", "element"), sep = "_", remove = TRUE)
  
  # TIDY HERBIVORES
  res.herbivores <- tidy(mode.herbivores, conf.int = TRUE, conf.level = 0.95) %>%
    filter(term != "(Intercept)") %>%
    mutate(full_name = elm) %>%
    separate(full_name, into = c("medium", "element"), sep = "_", remove = TRUE)
  
  # COMBINE
  results_list.global[[k]] <- bind_rows(res.water, res.herbivores)
  k <- k + 1
}

model_results.global <- bind_rows(results_list.global)

# --- DATA PREP ---
plot_data_global <- model_results.global %>%
  mutate(
    # Create Significance Column
    p.value = round(p.value, 3),
    sig_alpha = ifelse(p.value < 0.05, "Significant", "NS"),
    
    # Clean Labels for the Plot
    Treatment_Label = case_when(
      term == "Water.treatmentReduced" ~ "Reduced Amount",
      term == "Water.treatmentAltered Frequency" ~ "Reduced Frequency",
      term == "Herbivores1" ~ "Herbivores Removed", 
      TRUE ~ term
    ),
    
    # Create Panels (Water vs Herbivory)
    Panel = case_when(
      str_detect(term, "Water") ~ "Water Treatment",
      str_detect(term, "Herbivore") ~ "Herbivory",
      TRUE ~ "Other"
    ),
    
    # Set Factors for Ordering (Soil on top, Water on top)
    medium = factor(medium, levels = c("Soil", "Plant")),
    Panel = factor(Panel, levels = c("Water Treatment", "Herbivory"))) %>%
  mutate(element = factor(element, levels = sort(unique(element), decreasing = TRUE)))

# --- PLOT ---
ggplot(plot_data_global, aes(x = estimate, y = element, color = Treatment_Label, group = Treatment_Label)) +
  
  # The Reference Line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # The Points & Error Bars
  geom_pointrange(
    aes(xmin = conf.low, xmax = conf.high, alpha = sig_alpha),
    size = 0.3,
    position = position_dodge(width = 0.6),
    flatten = 0.1) +
  
    facet_grid(Panel + medium ~ ., scales = "free_y") +
  
  # Colors & Styling
  scale_color_manual(values = c(
    "Reduced Amount"    = "#D55E00",
    "Reduced Frequency" = "#0072B2",
    "Herbivores Removed"= "#009E73")) +
  
  scale_alpha_manual(values = c("NS" = 0.3, "Significant" = 1.0)) +
  
  labs(
    title = "",
    x = "Standardized Effect Size (ES)",
    y = NULL,
    color = "Treatment",
    alpha = "Significance"
  ) +
  
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    axis.title = element_text(size = 8, color = "black"),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text = element_text(color = "black", size = 6),
    legend.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 6, color = "black", face = "bold")
  )

ggsave("exports/DRIGrass_Fig.SUP_Global.tiff", width = 15, height = 15, units = "cm", dpi = 1200)

