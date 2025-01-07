##This file contains code to:
# - create multi-year df
# - prep Posthumas data frame and extract relevant SSDs
# - convert concentrations to TUs and create a sum_TU column
# - select chemicals that co-occur over multiple years
install.packages("pwr")
install.packages("reshape2")
install.packages("viridis")
install.packages("conflicted")

library(pwr)
library(reshape2)
library(viridis)
library(igraph)
library(dplyr)
library(data.table) 
library(conflicted)
library(tidyverse)
library(igraph)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

CONC_I2M2 <- DF

#Get Posthuma data ready for extracting SSD for toxic unit ----
POST_clean <- prio_posthuma_tox %>%
  select(CAS, Substance, `10LogSSDMedianConcentration(ug/L)-MuAcute EC50`, 
         `10LogSSDSlope(ug/L)-SigmaAcute EC50`, 
         `Remark on acute SSD slope`,
         `10LogSSDMedianConcentration(ug/L)-MuChronic NOEC`,
         `10LogSSDSlope(ug/L)-SigmaChronic NOEC`,
         `Remark on chronic SSD slope`)

unique(CHEMDF_w_max$CAS_Number)

POST_filtered <- POST_clean %>%
  filter(CAS %in% CHEMDF_w_max$CAS_Number)

##renaming for convenience
colnames(POST_filtered)[colnames(POST_filtered) == '10LogSSDMedianConcentration(ug/L)-MuAcute EC50'] <- 
  'LogSSDMedianConcentration_MuAcute_EC50'
colnames(POST_filtered)[colnames(POST_filtered) == '10LogSSDSlope(ug/L)-SigmaAcute EC50'] <- 
  'LogSSDSlope_SigmaAcute_EC50'
colnames(POST_filtered)[colnames(POST_filtered) == '10LogSSDMedianConcentration(ug/L)-MuChronic NOEC'] <- 
  'LogSSDMedianConcentration_MuChronic_NOEC'
colnames(POST_filtered)[colnames(POST_filtered) == '10LogSSDSlope(ug/L)-SigmaChronic NOEC'] <- 
  'LogSSDSlope_SigmaChronic_NOEC'

colnames(POST_filtered) <- gsub(" ", "_", colnames(POST_filtered), fixed = TRUE)

##use remarks to edit LogSSDSlope Posthuma
POST_filtered$LogSSDSlope_SigmaAcute_EC50[is.na(POST_filtered$LogSSDSlope_SigmaAcute_EC50) 
                                          & grepl("0.7", POST_filtered$Remark_on_acute_SSD_slope)] <- 0.7
POST_filtered$LogSSDSlope_SigmaChronic_NOEC[is.na(POST_filtered$LogSSDSlope_SigmaChronic_NOEC) 
                                            & grepl("0.7", POST_filtered$Remark_on_chronic_SSD_slope)] <- 0.7

## add non-log-tranformed SSD values to POST_filtered
POST_filtered <- POST_filtered %>%
  mutate(SSDMedianConcentration_MuChronic_NOEC = 10^LogSSDMedianConcentration_MuChronic_NOEC)
POST_filtered <- POST_filtered %>%
  mutate(SSDMedianConcentration_MuAcute_EC50 = 10^LogSSDMedianConcentration_MuAcute_EC50)

ssd_chronic <- setNames(POST_filtered$SSDMedianConcentration_MuChronic_NOEC, POST_filtered$CAS)
ssd_acute <- setNames(POST_filtered$SSDMedianConcentration_MuAcute_EC50, POST_filtered$CAS)
ssd_chronic_log<- setNames(POST_filtered$LogSSDMedianConcentration_MuChronic_NOEC, POST_filtered$CAS)
ssd_acute_log <- setNames(POST_filtered$LogSSDMedianConcentration_MuAcute_EC50, POST_filtered$CAS)


# Convert concentration values in chemical data to toxic units -----
#prep
colnames(CONC_I2M2)[-1] <- gsub("\\s+", "", colnames(CONC_I2M2)[-1])  
colnames(CONC_I2M2)[colnames(CONC_I2M2) == 'ResIndiceResultatBiologique'] <- 'I2M2'
colnames(CONC_I2M2)[colnames(CONC_I2M2) == 'year'] <- 'Year'
colnames(CONC_I2M2)[colnames(CONC_I2M2) == 'Sample_Source_Name.x'] <- 'Sample_Source_Name'

CONC_I2M2_clean <- CONC_I2M2 %>% #check which cols to select when doing real analysis
  select(Sample_Source_Name:`50-28-2`, I2M2)

long_CONC_I2M2 <- CONC_I2M2_clean %>%
  pivot_longer(
    cols = -c(Sample_Source_Name, Sample_Source_Code, Year, I2M2),
    names_to = "Chemical",
    values_to = "Concentration"
  )

# Calculate TUs and add as a new column
long_CONC_TU_I2M2 <- long_CONC_I2M2 %>%
  mutate(TU = map2_dbl(Chemical, Concentration, ~ {
    cas <- .x  
    concentration <- .y 
    chronic <- ssd_chronic[cas]
    acute <- ssd_acute[cas]

    if (!is.na(chronic)) {
      concentration / chronic
    } else if (!is.na(acute)) {
      concentration / acute
    } else {
      NA_real_ 
    }
  })) %>%
  filter(TU > 0)

#Analyse correlations between chemical concentrations 
head(long_CONC_TU_I2M2)
head(POST_filtered)

#aggregated across years
colnames(CONC_I2M2_clean) <- trimws(colnames(CONC_I2M2_clean))
POST_filtered$CAS <- trimws(POST_filtered$CAS)
cas_to_substance <- setNames(POST_filtered$Substance, POST_filtered$CAS)

CONC_I2M2_named <- CONC_I2M2_clean

colnames(CONC_I2M2_named)[colnames(CONC_I2M2_named) %in% names(cas_to_substance)] <- 
  cas_to_substance[colnames(CONC_I2M2_named)[colnames(CONC_I2M2_named) %in% names(cas_to_substance)]]

## correlation matrix all years
cor_matrix <- cor(CONC_I2M2_named %>% select(-c(Sample_Source_Name, Sample_Source_Code, Year, I2M2)),
                  use = "pairwise.complete.obs")

## visualise
cor_melted <- melt(cor_matrix)
chemical_order <- colnames(cor_matrix)

# Ensure that the 'Var1' and 'Var2' factors in cor_melted follow the same order as the correlation matrix
cor_melted$Var1 <- factor(cor_melted$Var1, levels = chemical_order)
cor_melted$Var2 <- factor(cor_melted$Var2, levels = chemical_order)

colnames(cor_melted) <- c("Var1", "Var2", "Correlation")

cor_heatplot <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  coord_fixed() +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 1.5) +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick", midpoint = 0) +  
  labs(title = "Correlation matrix of chemical concentrations", 
       x = "Chemical", 
       y = "Chemical") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8, face = "italic"), 
    axis.text.y = element_text(angle = 0, size = 7)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))  

#by year
numeric_data <- CONC_I2M2_named %>%
  mutate(across(-c(Sample_Source_Name, Sample_Source_Code, I2M2, Year), as.numeric))

cor_matrices <- numeric_data %>%
  select(-Sample_Source_Name, -Sample_Source_Code, -I2M2) %>%  
  group_by(Year) %>%  
  nest() %>%  
  mutate(
    cor_matrix = map(data, ~cor(.x %>% select(where(is.numeric)), use = "pairwise.complete.obs")) 
  ) %>%
  select(Year, cor_matrix)  

cor_long <- cor_matrices %>%
  mutate(cor_matrix_long = map(cor_matrix, ~as.data.frame(as.table(.)))) %>%
  select(Year, cor_matrix_long) %>%
  unnest(cor_matrix_long)

colnames(cor_long) <- c("Year", "Var1", "Var2", "Correlation")
cor_long$Var1 <- factor(cor_long$Var1, levels = chemical_order)
cor_long$Var2 <- factor(cor_long$Var2, levels = chemical_order)

cor_year_heatplot <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick", midpoint = 0) +
  facet_wrap(~Year, scales = "fixed") +  
  labs(title = "Correlations between chemicals per year",
       x = NULL,  
       y = NULL) + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())

print(cor_year_heatplot)

# attempt at network analysis ----

#presence-absence data site-year 
head(CONC_I2M2_named)

presence_absence <- CONC_I2M2_named %>%
  mutate(Site_Year = paste(Sample_Source_Code, Year, sep = "-")) %>%
  pivot_longer(cols = -c(Site_Year, Sample_Source_Name, Sample_Source_Code, Year, I2M2), 
               names_to = "Chemical", values_to = "Concentration") %>%
  mutate(Presence = ifelse(!is.na(Concentration) & Concentration > 0, 1, 0)) %>%
  select(Site_Year, Chemical, Presence) %>%
  pivot_wider(names_from = Chemical, values_from = Presence)

#matrix 
matrix <- presence_absence %>%
  select(-Site_Year) %>%  
  as.matrix()  

co_occurrence_matrix <- t(matrix) %*% matrix #matrix multipl. transpose
co_occurrence_matrix

diag(co_occurrence_matrix) <- 0 #correlation with itself
co_occurrence_matrix[is.na(co_occurrence_matrix)] <- 0

#create network and plot
graph <- graph_from_adjacency_matrix(co_occurrence_matrix, mode = "undirected", weighted = TRUE)

#distinguish communities - what makes most sense
community <- cluster_louvain(graph)
community <- cluster_edge_betweenness(graph)
community <- cluster_walktrap(graph)

plot(graph, 
     vertex.size = degree(graph),     # Node size based on degree
     edge.width = E(graph)$weight / 2000, # Scale edge width for visibility
     vertex.color = membership(community), # Color nodes by community
     layout = layout_with_fr(graph), # Circular layout
     main = "Chemical Co-occurrence Network with Communities")

# Assuming 'graph' is your graph object with weights
graph_layout <- layout_with_fr(graph, weights = E(graph)$weight)

community_colors <- c("tomato1", "plum3")
node_colors <- community_colors[membership(community)]

# Plot the graph
plot(graph, 
     layout = graph_layout,   
     vertex.size = degree(graph),
     vertex.color = node_colors,
     edge.width = E(graph)$weight / 2000,
     edge.color = "darkgray",  
     vertex.label = V(graph)$name,           # Use node names as labels
     vertex.label.color = "black",           # Change label color
     vertex.label.cex = 0.8, 
     main = "Network of chemicals")


summary(graph)
E(graph)$weight
# Calculate co-occurrence matrix
co_occurrence <- t(chemical_presence) %*% chemical_presence  # Transpose x presence/absence
diag(co_occurrence) <- 0  # Remove self-loops


#plot occurences ----
#plot chemical occurence over the years
dist_year <- long_TU_I2M2 %>%
  group_by(Sample_Source_Code, Year) %>%  
  summarise(chemical_count = n(), .groups = "drop") %>% 
  group_by(Year, chemical_count) %>%            
  summarise(sites = n(), .groups = "drop")    
dist <- long_TU_I2M2 %>%
  group_by(Sample_Source_Code) %>%              
  summarise(chemical_count = n(), .groups = "drop") %>%
  summarise(sites = n(), .groups = "drop") 

dist_plot_year <- ggplot(dist_year, aes(x = chemical_count, y = sites, color = as.factor(Year))) +
  geom_line() +
  geom_point() +
  labs(title = "Distribution per year",
       x = "Number of Chemicals", y = "Number of Sites", color = "Year") +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal()

dist_plot_year_violin <- ggplot(dist_year, aes(x = as.factor(Year), y = chemical_count, fill = as.factor(Year))) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Distribution of Chemicals per Site over Years",
       x = "Year", y = "Number of Chemicals", fill = "Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dist_plot_year_violin

dist_plot <- ggplot(dist, aes(x = chemical_count, y = sites)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Chemicals", y = "Number of Sites") +
  theme_minimal()

print(dist_plot_year)
print(dist_plot) 


# sum TU of selected chemicals----
selected_chems <- as.character(combo_count$Chemical_Combination[which.max(combo_count$n)])

selected_chems <- unlist(strsplit(selected_chems, ","))
selected_chems <- gsub('[\"()]', '', selected_chems)
selected_chems[1] <- sub("^c", "", selected_chems[1])
selected_chems <- trimws(selected_chems)

TU_I2M2_filtered <- TU_I2M2 %>%
  filter(
    rowSums(!is.na(TU_I2M2[selected_chems]) & TU_I2M2[selected_chems] != 0) == length(selected_chems)
  ) %>%
  select(Sample_Source_Name, Year, all_of(selected_chems))

#add column with sum of TUs
TU_I2M2 <- CONC_I2M2 %>% #wide format TU
  mutate(across(.cols = -c(Sample_Source_Name, Sample_Source_Code, I2M2, Year),
                .fns = ~ {
                  cas <- cur_column() 
                  chronic_value <- ssd_chronic[cas]
                  acute_value <- ssd_acute[cas]
                  
                  if (is.na(chronic_value)) {
                    if (!is.na(acute_value)) {
                      . / acute_value
                    } else {
                      NA
                    }
                  } else {
                    . / chronic_value
                  }
                }))


TU_I2M2 <- TU_I2M2 %>%
  dplyr::select(Sample_Source_Name:I2M2)

TU_I2M2 <- TU_I2M2 %>%
  select(-`50-28-2`)

TU_I2M2_sum <- TU_I2M2 %>% 
  mutate(Sum_TU = rowSums(select(., -Sample_Source_Name, -Sample_Source_Code, -Year, -I2M2), na.rm = TRUE))


#get long TU df
long_TU_I2M2 <- TU_I2M2_sum %>%
  select(Year, Sample_Source_Code, Sum_TU, I2M2) %>% 
  filter(Sum_TU > 0) 
 
#add variables from other dfs (environemental etc.)
# <- merge(DF, another_df[, c("Sample_Source_Name", "columnname")], by = "Sample_Source_Name", all.x = TRUE)
DF1 <- long_TU_I2M2
#add id 
DF1$id <- seq_along(DF1$Sample_Source_Code)

DF1 <- DF1[, c("id", setdiff(names(DF1), "id"))]

fwrite(CONC_I2M2, "CONC_I2M2.csv")
fwrite(TU_I2M2_sum, "TU_I2M2_sum.csv")
fwrite(long_CONC_TU_I2M2, "long_CONC_TU_I2M2.csv")
fwrite(long_TU_I2M2, "sum_TU.csv")
