library(dplyr)
library(tibble)
library(purrr)
library(adespatial)
library(tidyr)
library(ggplot2)
library(readr) 
library(cowplot)
library(vegan)

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

### look at cardoso ##

cardoso_list <- sp_degree_df %>% 
  filter(site_code == 'CARDC') 

cardoso_mat <- do.call(rbind, lapply(1:nrow(cardoso_list), FUN = function(x) data.frame(bromeliad_id = rep(cardoso_list$bromeliad_id[x], 
                                                                                     length(cardoso_list$composition[x])), 
                                                                  sp_id = unlist(cardoso_list$composition[x])))) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = "sp_id", values_from = 'presence', values_fill = 0) %>% 
  column_to_rownames("bromeliad_id") 

cardoso_beta <- beta.div(cardoso_mat, method = 'hellinger')

species_contribution_beta <- cardoso_beta$SCBD

subset_char <- function(namelist, char_vec){
  # check that all elements of namelist are in dimnames of matrix
  aux <- match(namelist, names(char_vec))
  sub_mat <- char_vec[aux]
  return(sub_mat)
}

species_contribution_beta_df <- cardoso_list %>% 
  select(bromeliad_id) %>% 
  rowwise() %>% 
  mutate(species_contribution_beta = list(species_contribution_beta)) 

cardoso_sp_spbd <- cardoso_list %>% 
  left_join(species_contribution_beta_df) %>% 
  mutate(spbd = map2(composition, species_contribution_beta, subset_char))

cardoso_sp_spbd$sp_degree[1]
cardoso_sp_spbd$spbd[1]

data.frame(sp_degree = cardoso_sp_spbd$sp_degree[[1]], spbd = cardoso_sp_spbd$spbd[[1]]) %>% 
  ggplot(aes(x = sp_degree, y = spbd)) + geom_point()


###### site level degree ####


site_level_degree <- lapply(prob_mat_list, colSums)



## beta diversity for all sites 

site_split <- split(sp_degree_df, sp_degree_df$site_code)

calculate_scbd <- function(df){
  
  #### calculates beta diversity based on presence absence, does not take into account abundance
  

  site_mat <- do.call(rbind, lapply(1:nrow(df), FUN = function(x) data.frame(bromeliad_id = rep(df$bromeliad_id[x], 
                                                                                                             length(df$composition[x])), 
                                                                                          sp_id = unlist(df$composition[x])))) %>% 
    mutate(presence = 1) %>% 
    pivot_wider(names_from = "sp_id", values_from = 'presence', values_fill = 0) %>% 
    column_to_rownames("bromeliad_id") 
  
  site_beta <- beta.div(site_mat, method = 'hellinger')
  
  species_contribution_beta <- site_beta$SCBD
  
  return(species_contribution_beta)
  
}

all_site_scbd <- lapply(site_split, calculate_scbd)



scbd_df <- all_site_scbd %>% 
  map(~data.frame(scbd = .x)) %>% 
  map_df(~ rownames_to_column(.x, "species_name"), .id= 'site')


degree_df <- site_level_degree %>% 
  map(~data.frame(degree = .x)) %>% 
  map_df(~ rownames_to_column(.x, "species_name"), .id= 'site')


scbd_df %>% 
  left_join(degree_df) %>% 
  ggplot(aes(x = degree, y = scbd)) +
  facet_wrap(~site) +
  geom_point() +
  geom_smooth(method = 'lm')


degree_df %>% 
  ggplot(aes(x = degree, fill = site)) +
  geom_histogram() +
  facet_wrap(~site)  



##### look at environmental gradients #######

trimmedbrom_pc <- read_csv(file = "../empirical_food_webs/05_models/data/trimmedbrom_pc.csv") %>% 
  mutate(visit_id = as.factor(visit_id),
         bromeliad_id = as.character(bromeliad_id),
         speciesPool = as.factor(speciesPool)) %>% 
         #site_code = as.factor(site_code)) %>% 
  filter(!is.na(max_water_combined)) %>% 
  filter(!is.na(total_detritus)) %>% 
  filter(!is.na(open.canopy))


max_water_all <- trimmedbrom_pc %>% 
  select(max_water_combined) 
  
trimmedbrom_pc %>% 
  ggplot(aes(x = max_water_combined)) +
  facet_wrap(~site_code) +
  geom_histogram() +
  scale_x_log10() +
  geom_histogram(data = max_water_all, aes(x = max_water_combined), alpha = 0.5) +
  theme_cowplot()
  
####### do mantel test ######

## environment dist #
bromleiads_environment <- trimmedbrom_pc %>% 
  select(bromeliad_id, site_code, log_detritus, max_water_combined) %>% 
  mutate(log_water = log(max_water_combined)) %>% 
  select(-max_water_combined)

bromeliads_split <- split(bromleiads_environment, bromleiads_environment$site_code)

bromeliad_enviornment_split <- bromeliads_split %>% 
  map(~select(.x, -site_code)) %>% 
  map(~column_to_rownames(.x, "bromeliad_id"))


bromeliad_environment_dist <- bromeliad_enviornment_split %>% 
  map(~dist(.x))

## beta diversity ##

abundance_correct_bromeliad <- abd_data %>% 
  filter(bromeliad_id %in% bromleiads_environment$bromeliad_id) %>% 
  left_join(site_codes) %>% 
  select(bwg_name, bromeliad_id, abundance, site_code)

abundance_split <- split(abundance_correct_bromeliad, abundance_correct_bromeliad$site_code)

abundance_wider <- abundance_split %>% 
  map(~select(.x, -site_code)) %>% 
  map(~group_by(.x, bwg_name, bromeliad_id)) %>% 
  map(~summarise(.x, abundance_total = sum(abundance))) %>% 
  map(~pivot_wider(.x, names_from = bwg_name, values_from = abundance_total, values_fill = 0)) %>% 
  map(~column_to_rownames(.x, "bromeliad_id"))

beta_div_abundance <- abundance_wider %>% 
  map(~vegdist(.x, binary=FALSE, method = "bray"))


### mantel test ##

mantel_loo_function <- function(site) {
  
  mantel_res <- mantel(beta_div_abundance[[site]], bromeliad_environment_dist[[site]])

  mantel_loo_all <- c()
  
  for(i in 1:ncol(abundance_wider[[site]])){
    
    new_abundance <- abundance_wider[[site]][, -i]
    
    beta_loo <- vegdist(new_abundance, binary=FALSE, method = "bray")
    
    mantel_loo <- mantel(beta_loo, bromeliad_environment_dist[[site]])
    
    mantel_loo_all <- c(mantel_loo_all, mantel_loo$statistic)
    
  }
  
  mantel_all_df <- data.frame(bwg_name =   colnames(abundance_wider[[site]]), mantel_r_loo = mantel_loo_all, mantel_original_r = mantel_res$statistic, site_code = site)
  
  return(mantel_all_df)
}

mantel_loo_function(names(beta_div_abundance)[2])

mantel_all_sites <- lapply(names(beta_div_abundance), safely(mantel_loo_function))
## need to figure out why 4 sites didn't work. safely is great!

mantel_all_sites_df_all <- mantel_all_sites %>% 
  map(~ .x$result) %>% 
  keep(~ is.null(.x) == FALSE) %>% 
  map_df(~data.frame(.x)) %>% 
  mutate(mantel_ratio = mantel_r_loo/mantel_original_r)

site_degree_df <- site_level_degree %>% 
  map(~data.frame(degree = .x)) %>% 
  map_df(~rownames_to_column(.x, "bwg_name"),  .id = 'site_code')

mantel_contribution_degree <- mantel_all_sites_df_all %>% 
  left_join(site_degree_df) %>% 
  group_by(site_code) %>% 
  mutate(sd_mantel_loo = sd(mantel_r_loo)) %>%
  mutate(z_score = (mantel_original_r - mantel_r_loo)/sd_mantel_loo)


mantel_contribution_degree %>% 
  ggplot(aes(x = degree, y = z_score)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~site_code, scales = 'free') 

mantel_contribution_degree %>% 
  ggplot(aes(x = mantel_r_loo)) +
  geom_histogram() + 
  geom_vline(aes(xintercept = mantel_original_r), colour = 'red') +
  facet_wrap(~site_code, scales = 'free') 

library(ggridges)
mantel_contribution_degree %>% 
  ungroup %>% 
  mutate(site_code_f = forcats::fct_reorder(site_code, mantel_r_loo)) %>% 
  ggplot(aes(x = mantel_r_loo, y = site_code_f)) + 
  stat_density_ridges()
  
water_site <- trimmedbrom_pc %>% 
  group_by(site_code) %>% 
  summarise(sd_wat = sd(max_water_combined), n = n())

mantel_contribution_degree %>% 
  left_join(water_site) %>% 
  ggplot(aes(x = sd_wat, y = mantel_r_loo)) +
  geom_point() +
  geom_point(aes(y = mantel_original_r), size = 2, col = "green") + 
  geom_hline(yintercept = 0)

water_site %>% 
  filter(sd_wat > 1e3)

site_order <- water_site %>% 
  arrange(sd_wat) %>% select(site_code) %>% unlist()

trimmedbrom_pc %>% 
  filter(site_code == "CABO") %>% View()




mantel_contribution_degree %>% 
  mutate(site_code = factor(site_code, levels = site_order)) %>% 
  ggplot(aes(x = z_score, y = degree)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~site_code, scales = 'free') 

mantel_contribution_degree %>% 
  mutate(site_code = factor(site_code, levels = site_order)) %>% 
  ggplot(aes(x = mantel_r_loo)) +
  geom_histogram() + 
  geom_vline(aes(xintercept = mantel_original_r), colour = 'red') +
  facet_wrap(~site_code, scales = 'free') 




#### get better food web metrics###

library(cheddar)
library(igraph)


prob_mat_list[[1]]






1#### check errors in data ###
abundance_wider %>% 
  keep(~any(map_lgl(.x, ~typeof(.x) == "list"))) 

#### Cardoso ###

# abundance_wider$CARDO
# 
# table(abundance_split$CARDO$bwg_name, abundance_split$CARDO$bromeliad_id)
#   
# pivot_wider(names_from = bwg_name, values_from = abundance)
# 
# abd_data %>% 
#   filter(location == "Cardoso_open") %>% 
#   filter(bwg_name == "Coleoptera.48") %>% View()
# 
# 
# abd_data %>% 
#   filter(bwg_name == "Coleoptera.48") %>% 
#   count(species_id, bromeliad_id, location) %>% 
#   View()
# 
# 
# 
# #### EVTA ###
# 
# abundance_wider$EVTA 
# 
# abundance_split$EVTA %>% 
#   count(bwg_name, bromeliad_id) %>% 
#   filter(n>1)
# 
# table(abundance_split$EVTA$bwg_name, abundance_split$EVTA$bromeliad_id)
# 
# pivot_wider(names_from = bwg_name, values_from = abundance)
# 
# abd_data %>% 
#   filter(location == "Elverde_tabunoco") %>% 
#   filter(bwg_name == "Coleoptera.19") %>% View()
# 
# abd_data %>% 
#   filter(bwg_name == "Coleoptera.19") %>% 
#   count(species_id, bromeliad_id, location) %>% 
#   View()
