# Data cleaning and organisation for DRI-GRASS publication:
# Clean environment:
rm(list=ls())

# Load packages:
library(tidyverse)
library(httr2)

# Data Import (via Figshare API) --------------------------------------------

# Define Figshare article ID containing the pigment dataset
article_id <- 29516558

# Request metadata for the specified Figshare article
res <- request(paste0("https://api.figshare.com/v2/articles/", article_id)) %>%
  req_perform()

# Parse the JSON response into an R list / data structure
meta <- res %>%
  resp_body_json(simplifyVector = TRUE)

# Select the specific file (first file associated with the record).
file_url <- meta$files$download_url[1] 
dat <- read_csv(file_url, show_col_types = FALSE)

# Cleaning:
datos <- dat %>%
  mutate(Year = as.factor(Year),
         Water.treatment = as.factor(Water.treatment),
         Herbivores = as.factor(Herbivores),
         Plot.number = as.factor(Plot.number)) %>%
  select(-c(Open.control)) %>%
  mutate(WTH = paste(Water.treatment, Herbivores, sep = ".")) %>%
  mutate(WHS = paste(Water.treatment, Herbivores, Year, sep = ".")) %>%
  mutate(WTH = recode_factor(WTH,
                             "Control.0" = "Ambient -RH",
                             "Control.1" = "Ambient +RH",
                             "Reduced.0" = "Reduced amount -RH",
                             "Reduced.1" = "Reduced amount +RH",
                             "Altered Frequency.0" = "Reduced frequency -RH",
                             "Altered Frequency.1" = "Reduced frequency +RH"),
         WTH = factor(WTH, levels = c(
           "Ambient -RH", "Ambient +RH", "Reduced amount -RH", "Reduced amount +RH",
           "Reduced frequency -RH", "Reduced frequency +RH"))) %>%
  mutate(Plot.number = as.factor(Plot.number),
         Treatment = as.factor(Treatment),
         Water.treatment = as.factor(Water.treatment),
         Herbivores = as.factor(Herbivores),
         WHS = as.factor(WHS),
         WTH = as.factor(WTH),
         Year = as.factor(Year)) %>%
  select(Plot.number, Treatment, Water.treatment, WHS, Herbivores, WTH, Year, everything())

# Export:
write.csv(datos, "data/DRIGrass_tidy.data.csv", row.names = FALSE)