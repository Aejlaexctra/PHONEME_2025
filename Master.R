# --- SET RANDOM SEED ---
set.seed(42)

# --- (1) DOWNLOADING DATASETS AND MODULES --- ####

## --- Import necessary packages --- ####
library(nloptr)
library(nlme)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)

## --- Loading Phylogenetic, Spatial Distance and Contact Matrices --- ####
phylomatrix <- read.csv("data/distance_matrices/phylogenetic_covariance_matrix.csv",
                        check.names = FALSE)
rownames(phylomatrix) <- phylomatrix[,1]
phylomatrix <- phylomatrix[,-1]
phylomatrix <- as.matrix(phylomatrix)
spmatrix <- read.csv("data/distance_matrices/spatial_distance_matrix.csv",
                     check.names = FALSE)
rownames(spmatrix) <- spmatrix[,1]
spmatrix <- spmatrix[,-1]
spmatrix <- as.matrix(spmatrix)

## --- Loading NEE22 (3) --- ####
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)
colnames(nee22)[colnames(nee22) == "island"] <- "Island.Endemic"

## --- Loading NEE24 (1) --- ####
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)
colnames(nee24)[colnames(nee24) == "Island.Endemic"] <- "Island.Endemic"

## --- Loading Phoible (4) --- ####
var_phon_inv_spec_data <- read.csv("Data/cldf-datasets-inventory-study/PHOIBLE-data.csv")
var_phon_inv_glot_data <- read.csv("Data/cldf-datasets-inventory-study/phoible/cldf/languages.csv")
phoible <- merge(
  var_phon_inv_spec_data[, c("Glottocode", "Sounds", "Latitude", "Longitude")],
  var_phon_inv_glot_data[, c("ISO639P3code", "Glottocode")], by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)
colnames(phoible)[colnames(phoible) == "Sounds"] <- "Phoneme.Inventory.Size"

# --- (2) FUNCTION DEFINITIONS --- ####

## --- GLS (2) Setup --- ####
best_p <- function (p,formula,data,spmatrix,phylomatrix) {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  if (inherits(res,"try-error")) {
    out <- -10000
  } else {
    out <- res$logLik
    # if (res$logLik>0) {out <- -10000}
  }
  -out
}

ml_fit <- function (p,formula,data,spmatrix,phylomatrix) {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  res
}

# --- (3) CLEANING AND FORMATING DATASETS --- ####

## --- A1a --- ####

# --- Merging Datasets A1a ---
A1a <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region", "Island.Endemic")],
  phoible[, c("ISO639P3code", "Glottocode", "Latitude", "Longitude", "Phoneme.Inventory.Size")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)

# --- Data Cleanup A1a ---
# Remove all NA entries
print("Size of initial dataset:")
print(dim(A1a)[1])
old_A1a <- A1a
A1a <- na.omit(A1a) 
print("Size of dataset with NA removed:")
print(dim(A1a)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- which(table(A1a$ISO) != 1)
print(which(table(A1a$ISO) != 1))
A1a <- A1a[!(A1a$ISO %in% names(which(table(A1a$ISO) != 1))),]
print("Size of dataset after duplicates removed:")
print(dim(A1a[1]))

# --- Adjusting data with matrices A1a ---
common_ids <- Reduce(intersect, list(
  A1a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A1a <- A1a[A1a$ISO %in% common_ids, ]
A1a_phylomatrix <- phylomatrix[common_ids, common_ids]
A1a_spmatrix <- spmatrix[common_ids, common_ids]

## --- A1b --- ####

# --- Merging Datasets A1b ---
A1b <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "Island.Endemic")], by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Data Cleanup A1b ---
# Remove all NA entries
print("Size of initial dataset:")
print(dim(A1b)[1])
old_A1a <- A1b
A1b <- na.omit(A1b) 
print("Size of dataset with NA removed:")
print(dim(A1b)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- which(table(A1b$ISO) != 1)
print(which(table(A1b$ISO) != 1))
A1b <- A1b[!(A1b$ISO %in% names(which(table(A1b$ISO) != 1))),]
print("Size of dataset after duplicates removed:")
print(dim(A1b[1]))

# --- Adjusting data with matrices A1b ---
common_ids <- Reduce(intersect, list(
  A1b$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A1b <- A1b[A1b$ISO %in% common_ids, ]
A1b_phylomatrix <- phylomatrix[common_ids, common_ids]
A1b_spmatrix <- spmatrix[common_ids, common_ids]

## --- A2 --- ####

## Merging Datasets
A2 <- merge(
  nee22[, 
        c("ISO", "region", "L1_pop", "bordering_language_richness", "altitude_range", "area")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "Island.Endemic")], 
  by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Data Cleanup ---
# Remove all NA entries
paste("Size of initial dataset:", dim(A2)[1])
old_A2 <- A2
A2 <- na.omit(A2) 
paste("Size of dataset with NA removed:", dim(A2)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- duplicated(A2$ISO)
print(data$ISO[iso_duplicates])
data <- A2[!unlist(iso_duplicates),]
paste("Size of dataset after duplicates removed:", dim(A2[1]))

# --- Adjusting data with matrices ---
common_ids <- Reduce(intersect, list(
  A2$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A2 <- A2[A2$ISO %in% common_ids, ]
A2_phylomatrix <- phylomatrix[common_ids, common_ids]
A2_spmatrix <- spmatrix[common_ids, common_ids]

## --- A3 --- ####

## --- A4 --- ####

# --- (4) EXPORT DATASETS FOR ANALYSIS --- ####

## --- A1a --- ####

# --- Preview and save data A1a ---
paste("Size of dataset adjusted:", dim(A1a)[1])
paste("#Island.Endemic: ", sum(A1a$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A1a$Island.Endemic) / dim(A1a)[1])
write.csv(A1a, file = "output/A1a/A1a_adjusted.csv", row.names=FALSE)

# --- Preview and save matrices A1a ---
print(dim(A1a_phylomatrix))
print(dim(A1a_spmatrix))
print(A1a_phylomatrix[1:10,1:10])
print(A1a_spmatrix[1:10,1:10])
save(A1a_phylomatrix, file = "output/A1a/A1a_phylomatrix.RData") 
save(A1a_spmatrix, file = "output/A1a/A1a_spmatrix.RData")

## --- A1b --- ####

# --- Preview and save data A1b ---
paste("Size of dataset adjusted:", dim(A1b)[1])
paste("#Island.Endemic: ", sum(A1b$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A1b$Island.Endemic) / dim(A1b)[1])
write.csv(A1b, file = "output/A1b/A1b_adjusted.csv", row.names=FALSE)

# --- Preview and save matrices A1b ---
print(dim(A1b_phylomatrix))
print(dim(A1b_spmatrix))
print(A1b_phylomatrix[1:10,1:10])
print(A1b_spmatrix[1:10,1:10])
save(A1b_phylomatrix, file = "output/A1b/A1b_phylomatrix.RData") 
save(A1b_spmatrix, file = "output/A1b/A1b_spmatrix.RData")

## --- A2 --- ####

# --- Preview and save data ---
paste("Size of dataset adjusted:", dim(A2)[1])
paste("#Island_Endemic: ", sum(A2$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(A2$Island.Endemic) / dim(A2)[1])
write.csv(A2, file = "output/A2/A2_adjusted.csv", row.names=FALSE)

### Preview and save matrices
print(dim(A2_phylomatrix))
print(dim(A2_spmatrix))
print(A2_phylomatrix[1:10,1:10])
print(A2_spmatrix[1:10,1:10])
save(A2_phylomatrix, file = "output/A2/A2_phylomatrix.RData") 
save(A2_spmatrix, file = "output/A2/A2_spmatrix.RData")

## --- A3 --- ####

## --- A4 --- ####

# --- (5) DATASET VISUALISATION --- ####

## Setup ####

region_order <- c("Oceania", "Australia and New Zealand", 
                  "South-Eastern Asia", "Southern Asia", "Asia", 
                  "Europe", "North Africa and Arabia", "Africa", 
                  "Western Africa", "Northern America", "Central America", "South America")

# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE) 
}

## --- A1a --- ####

# --- Get raw A1a adjusted data and matrices ---
raw_A1a <- read.csv("output/A1a/A1a_adjusted.csv")
load("output/A1a/A1a_phylomatrix.RData")
load("output/A1a/A1a_spmatrix.RData")

# --- Preview raw data A1a ---
head(raw_A1a)
summary(raw_A1a)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1a, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages Observed",
    x = "PISa",
  ) + 
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15)

ggsave(
  filename = "output/A1a/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Log transform variables A1a ---
A1a <- raw_A1a
A1a$Phoneme.Inventory.Size <- log(A1a$Phoneme.Inventory.Size) 
A1a$L1_pop <- log(A1a$L1_pop + 0.5)

# --- Mapping phoneme inventory sizes onto world map A1a ---
# Get world map data
world <- map_data("world")
PIS_world_map <-ggplot() +
  # Base world map
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray90", color = "gray50"
  ) +
  # Points with size and color based on count
  geom_point(
    data = A1a,
    aes(x = Longitude, y = Latitude, color = log(Phoneme.Inventory.Size)),
    alpha = 0.7
  ) +
  scale_color_viridis_c(option="plasma") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(x = "", y = "", color = "Log PISa")
ggsave("output/A1a/PIS_world_map.png", plot = PIS_world_map, width = 15, height = 10, dpi = 300)
PIS_world_map

# --- Preview correlations A1a ---
# Plot data
phoneme_L1_scatter = ggplot(A1a, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  # stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) + 
  labs(x = "L1_pop", y = "PISa", color = "Region") + # title = "Phoneme Count ~ L1 Population Size Scatterplot",  
  theme_minimal()
# theme(legend.position = "none") + # Hide the legend for now, often too many groups
# Save graph
ggsave("output/A1a/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias A1a ---
# number of repetitions
n_rep = 1e4
# Get region names
region_names <- region_order
# Find number of languages in each region for A1a
A1a_region_tally <- as.data.frame(table(A1a$region))
A1a_region_tally <- A1a_region_tally[match(region_names, A1a_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_names), nrow = n_rep))
colnames(region_sample_tallies) <- region_names # Set column names
# Repeated n_rep times, randomly sample without replacement 
for(i in 1:n_rep){
  # Get region sample
  region_sample <- sample(nee22$region, size = length(A1a$region), replace = FALSE)
  # Current region tally.
  region_sample_tally = table(region_sample)
  # Add current tally to all tally data
  region_sample_tallies[i, names(region_sample_tally)] <- as.numeric(region_sample_tally)
}
# Reshape the dataframe to long format
region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
region_sample_tallies_long$Group <- factor(region_sample_tallies_long$Group, levels = region_names)
region_box <- ggplot(region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1a_region_tally, 
             aes(x = region_names, y = Freq), 
             color = "red",
             size = 2) + 
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1a/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# --- Phoneme.Inventory.Size across regions A1a ---
box_phoneme_regions <- ggplot(raw_A1a, aes(
  x = factor(region, levels = region_order),
  y = Phoneme.Inventory.Size
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median(raw_A1a$Phoneme.Inventory.Size), color = "red") +
  labs(
    x = "Region",
    y = "PISa"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(box_phoneme_regions)
ggsave("output/A1a/phoneme_region_box.png", plot = box_phoneme_regions, width = 15, height = 15, dpi = 300, scale = 0.5)

## --- A1b --- ####

# --- Preview raw data A1b ---

# --- Get raw A1b adjusted data and matrices ---
raw_A1b <- read.csv("output/A1b/A1b_adjusted.csv")
load("output/A1b/A1b_phylomatrix.RData")
load("output/A1b/A1b_spmatrix.RData")

head(raw_A1b)
summary(raw_A1b)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1b, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages Observed",
    x = "PISb",
  ) + 
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15)

ggsave(
  filename = "output/A1b/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Log transform variables A1b ---
A1b<- raw_A1b
A1b$Phoneme.Inventory.Size <- log(A1b$Phoneme.Inventory.Size) 
A1b$L1_pop <- log(A1b$L1_pop + 0.5)

# --- Preview correlations A1b ---
# Plot data
phoneme_L1_scatter = ggplot(A1b, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  # stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) + 
  labs(x = "L1_pop", y = "PISb", color = "Region") + # title = "Phoneme Count ~ L1 Population Size Scatterplot",  
  theme_minimal()
# theme(legend.position = "none") + # Hide the legend for now, often too many groups
# Save graph
ggsave("output/A1b/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias A1b ---
# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE) 
}
# number of repetitions
n_rep = 1e4
# Get region names
region_names <- region_order
# Find number of languages in each region for A1b
A1b_region_tally <- as.data.frame(table(A1b$region))
A1b_region_tally <- A1b_region_tally[match(region_names, A1b_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_names), nrow = n_rep))
colnames(region_sample_tallies) <- region_names # Set column names
# Repeated n_rep times, randomly sample without replacement 
for(i in 1:n_rep){
  # Get region sample
  region_sample <- sample(nee22$region, size = length(A1b$region), replace = FALSE)
  # Current region tally.
  region_sample_tally = table(region_sample)
  # Add current tally to all tally data
  region_sample_tallies[i, names(region_sample_tally)] <- as.numeric(region_sample_tally)
}
# Reshape the dataframe to long format
region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
region_sample_tallies_long$Group <- factor(region_sample_tallies_long$Group, levels = region_names)
region_box <- ggplot(region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1b_region_tally, 
             aes(x = region_names, y = Freq), 
             color = "red",
             size = 2) + 
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1b/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# --- Phoneme.Inventory.Size across regions A1b ---
box_phoneme_regions <- ggplot(raw_A1b, aes(
  x = factor(region, levels = region_order),
  y = Phoneme.Inventory.Size
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median(raw_A1b$Phoneme.Inventory.Size), color = "red") +
  labs(
    x = "Region",
    y = "PISb"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(box_phoneme_regions)
ggsave("output/A1b/phoneme_region_box.png", plot = box_phoneme_regions, width = 15, height = 15, dpi = 300, scale = 0.5)

## --- A2 --- ####

# Get raw A2 adjusted data and matrices
raw_A2 <- read.csv("output/A2/A2_adjusted.csv")
load("output/A2/A2_phylomatrix.RData")
load("output/A2/A2_spmatrix.RData")

# --- Log transform variables ---
A2 <- raw_A2
A2$Phoneme.Inventory.Size <- log(A2$Phoneme.Inventory.Size) 
A2$area <- log(A2$area)
A2$altitude_range <- log(A2$altitude_range + 0.5)
A2$L1_pop <- log(A2$L1_pop + 0.5)
A2$bordering_language_richness <- log(A2$bordering_language_richness + 0.5)

# --- Preview transformed data ---
head(A2)
summary(A2)
PIS_Island <- ggplot(data = A2, aes(x = as.factor(Island.Endemic), y = Phoneme.Inventory.Size)) +
  geom_boxplot()
print(PIS_Island)

# --- Preview correlations ---
pairs(A2[,
         c("Phoneme.Inventory.Size", "area", "L1_pop", 
           "altitude_range", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A2[,
       c("Phoneme.Inventory.Size", "area", "L1_pop", 
         "altitude_range", "bordering_language_richness")])
cor.test(A2$Phoneme.Inventory.Size,A2$Island.Endemic)

## --- A3 --- ####

## --- A4 --- ####

# --- (6) OLS AND GLS ANALYSIS --- ####

## Setup ####
A1_model <- c(Phoneme.Inventory.Size ~ L1_pop)

## --- A1a --- ####

# --- Ordinary Least Squares A1a ---
A1a_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1a)
print(summary(A1a_l1))
print(BIC(A1a_l1))
print(logLik(A1a_l1))

# --- Best p values for all models A1a ---
A1a_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=A1_model[[1]],
                   data=A1a,
                   spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1a_p.res, file = "output/A1a/A1a_p_res.RData")

# --- Maximum likelihood fits for all models A1a ---
A1a_model <- ml_fit(p=A1a_p.res$par,formula=A1_model[[1]],data=A1a,spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix)
# Save model fit
save(A1a_model, file = "output/A1a/A1a_model.RData") 

## --- A1b --- ####

# --- Ordinary Least Squares A1b ---
A1b_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1b)
print(summary(A1b_l1))
print(BIC(A1b_l1))
print(logLik(A1b_l1))

# --- Best p values for all models A1b ---
A1b_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=A1_model[[1]],
                   data=A1b,
                   spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1b_p.res, file = "output/A1b/A1b_p_res.RData")

# --- Maximum likelihood fits for all models A1b ---
A1b_model <- ml_fit(p=A1b_p.res$par,formula=A1_model[[1]],data=A1b,spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix)
# Save model fit
save(A1b_model, file = "output/A1b/A1b_model.RData") 

## --- A2 --- ####

# --- Ordinary Least Squares ---
A2_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A2)
print(summary(A2_l1))
print(BIC(A2_l1))
print(logLik(A2_l1))

# --- Model predictor combinations ---
A2_predictors <- c("Island.Endemic", "area",
                "altitude_range", "bordering_language_richness",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A2_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A2_predictors, x))
A2_models <- c(Phoneme.Inventory.Size ~ 1)
for (i in 1:n_pred) {
  for (j in 1:dim(pred_combs[[i]])[2]) {
    combination <- pred_combs[[i]][,j]
    # Accounting for interactions
    # if(all(c("Island.Endemic","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "Island.Endemic:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    # if(all(c("L1_pop","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "L1_pop:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    # if(all(c("bordering_language_richness","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "bordering_language_richness:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    model <- (paste(response, paste(combination, collapse="+"), sep="~"))
    model <- as.formula(model)
    A2_models <- c(A2_models,model)
  }
}
unlist(A2_models)

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A2,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(p.res, file = "output/A2/ml_p_res.RData")
spmatrix_temp <- spmatrix/max(spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/p.res$par[2])^2)
mat <- as.matrix((p.res$par[3]*(1-p.res$par[1])*spmatrix_temp+(1-p.res$par[1])*(1-p.res$par[3])*phylomatrix+p.res$par[1]*diag(dim(phylomatrix)[1])))
save(mat, file = "output/A2/mat.RData")

# --- Maximum likelihood fits for all models ---
ml_fits <- mclapply(
  X = A2_models,
  FUN = function(f) ml_fit(p=p.res$par,formula=f,data=A2,spmatrix=spmatrix,phylomatrix=phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(ml_fits, "output/A2/ml_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(ml_fits)
ml_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "area" = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(ml_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A2)[1])
    next
  }
  # Get coefficient values
  coef_val <- ml_summary$tTable[,1]
  # Get corresponding p-values
  coef_p <- ml_summary$tTable[,4]
  # Get coefficient names
  coef_names <- names(ml_summary$tTable[,4])
  # Get number of coefficients
  n_coef <- length(coef_names)
  # Add coef info into model_summaries
  for (j in 2:n_coef){
    if(coef_p[j]<=0.05){
      ml_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
    } else {
      ml_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
    }
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A2)[1])
}

# Create new data frame from lists
ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(ml_summaries$Island.Endemic),
                      unlist(ml_summaries$area),
                      unlist(ml_summaries$altitude_range),
                      unlist(ml_summaries$bordering_language_richness),
                      unlist(ml_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(ml_data) <- c("BIC", "logLik", "MSE", A2_predictors)
# Sort based on increasing BIC
ml_data <- ml_data[order(ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(ml_data$BIC)
ml_data$deltaBIC <- ml_data$BIC - min_BIC
# Save output
write.csv(ml_data, file = "output/A2/ml_data.csv", row.names=FALSE)

## --- A3 --- ####

## --- A4 --- ####

# --- (7) REFERENCES --- ####
#
# (1) Bromham, L., Yaxley, K.J. & Cardillo, M. Islands are engines of language diversity. Nat Ecol Evol 8, 1991–2002 (2024). https://doi.org/10.1038/s41559-024-02488-4
# 
# (2) Hua, X., Greenhill, S.J., Cardillo, M. et al. The ecological drivers of variation in global language diversity. Nat Commun 10, 2047 (2019). https://doi.org/10.1038/s41467-019-09842-2
# 
# (3) Bromham, L., Dinnage, R., Skirgård, H. et al. Global predictors of language endangerment and the future of linguistic diversity. Nat Ecol Evol 6, 163–173 (2022). https://doi.org/10.1038/s41559-021-01604-y
# 
# (4) Cormac Anderson, Tiago Tresoldi, Simon J Greenhill, Robert Forkel, Russell Gray, Johann-Mattis List, Variation in phoneme inventories: quantifying the problem and improving comparability, Journal of Language Evolution, Volume 8, Issue 2, July 2023, Pages 149–168, https://doi.org/10.1093/jole/lzad011
# 