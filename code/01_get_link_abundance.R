library(dplyr)
library(tibble)
library(purrr)

## bromeliad data

site_codes <- 
  readr::read_csv("../empirical_food_webs/probabilistic_matrices/site_codes.csv")

# Get species composition of each bromeliad
abd_data <- 
  readr::read_csv("../empirical_food_webs/01_combined_food_web/data/abundance_cleaned.csv")

# bromeliad lookup table
brm_data <- 
  readr::read_csv("../empirical_food_webs/01_combined_food_web/data/bromeliads_cleaned.csv")
    

##### get webs ####

source("../empirical_food_webs/03_food_web_metrics/food_web_functions.R")

# Find the data folder using the `here` package
filenames <- 
  list_site_matrices2(data_folder = "../empirical_food_webs/03_food_web_metrics/data/pasted_matrices/")

# read in all probability matrices as a list
prob_mat_list <- read_prob_matrices(filenames)



site_bromeliad_level_composition <- readr::read_rds("../empirical_food_webs/03_food_web_metrics/data/site_bromeliad_level_composition.rds")


# store prob mat list in a data frame

prob_mat_df <- enframe(prob_mat_list,
                       name = "site_code",
                       value = "prob_mat")


sub_mats_df <- site_bromeliad_level_composition %>% 
  left_join(prob_mat_df) %>% 
  mutate(sub_mats = map2(composition, prob_mat, subset_mat))

##

sp_degree_df <- sub_mats_df %>% 
  mutate(sp_degree = map(sub_mats, colSums))

sp_degree_df %>% 
  filter(site_code == 'CARDC')
  
hist(unlist(sp_degree_df$sp_degree[1]))
