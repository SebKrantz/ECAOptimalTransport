#####################################################################
# Transport Network: Summary Statistics for Different Countries
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, units, tmap, ggplot2, install = TRUE)
fastverse_conflicts()

files <- list.files("results/grid_network/country", pattern = "_sigma38_") |> 
  grep(pattern = "_realloc_", invert = FALSE, value = TRUE) |> 
  grep(pattern = "_irs_", invert = FALSE, value = TRUE) |> 
  grep(pattern = "_ann_", value = TRUE, invert = TRUE) |>
  grep(pattern = "_ug_", value = TRUE, invert = TRUE)
shift <- 0
countries <- substr(files, 15+shift, 17+shift) |> unique()
setdiff(unique(substr(list.files("data/grid_network/country"), 1, 3)), countries)

iso3_codes <- c(
  "KGZ" = "Kyrgyzstan",
  "TJK" = "Tajikistan",
  "BGR" = "Bulgaria",
  "EST" = "Estonia",
  "UKR" = "Ukraine",
  "GRC" = "Greece",
  "POL" = "Poland",
  "SRB" = "Serbia",
  "TKM" = "Turkmenistan",
  "BLR" = "Belarus",
  "TUR" = "Turkey",
  "AZE" = "Azerbaijan",
  "UZB" = "Uzbekistan",
  "ROU" = "Romania",
  "ALB" = "Albania",
  "BIH" = "Bosnia",
  "ARM" = "Armenia",
  "MNE" = "Montenegro",
  "HRV" = "Croatia",
  "XKX" = "Kosovo"
)
setdiff(names(iso3_codes), countries)
files <- files[substr(files, 15+shift, 17+shift) %in% names(iso3_codes)]

results <- list(
  nodes = lapply(files[startsWith(files, "nodes")], function(f) fread(paste0("results/grid_network/country/", f))) |>
    set_names(substr(files[startsWith(files, "nodes")], 15+shift, 17+shift)),
  edges = lapply(files[startsWith(files, "edges")], function(f) fread(paste0("results/grid_network/country/", f))) |>
    set_names(substr(files[startsWith(files, "edges")], 15+shift, 17+shift))
)

summary_stats <- lapply(names(iso3_codes), function(c) {
  nodes <- results$nodes[[c]] |> subset(!is_buff)
  edges <- results$edges[[c]] |> subset(!is_buff)
  
  data.frame(
    country = iso3_codes[c],
    nodes = nrow(nodes),
    edges = nrow(edges),
    pop = fmedian(nodes$pop_cell),
    GDPC = fmedian(nodes$cell_GDPC_const_2017_PPP) * 1e6,
    TE = fmedian(edges$time_efficiency),
    delta_tau = fmedian(0.1158826 * log(edges$sp_distance)),
    delta_I = fmedian(edges$cost_kmh)
  )
}) |> rbindlist()

summary_stats |> 
  xtable::xtable(digits = 2) |> 
  print(include.r = FALSE, booktabs = TRUE)

summary_stats |> num_vars() |> pwcor(P = TRUE)


# Create Results Table for Paper: 
tab <- readxl::read_xlsx("results/grid_network/ALL_Mrealloc_table.xlsx", sheet = "EXP_IRS") |> 
       subset(ckmatch(names(iso3_codes), iso3c)) |> 
       join(readxl::read_xlsx("results/grid_network/ALL_Mrealloc_table.xlsx", sheet = "UG_IRS") |> 
            rename(exp_match_realloc_perc = ug_match_realloc_perc), on = "iso3c", drop = TRUE) |> 
       join(readxl::read_xlsx("results/grid_network/ALL_cent_edges_perc.xlsx", sheet = "Results") |> 
              select(iso3c = country, realloc_fixed_irs), on = "iso3c") 

tab |> 
  select(country, realloc_gain_perc, realloc_fixed_irs, exp_match_realloc_perc, ug_match_realloc_perc) |> 
  xtable::xtable(digits = 3) |> print(include.r = FALSE, booktabs = TRUE)
  

# Also add / recalculate cost of matching investments
costs <- lapply(names(iso3_codes), function(c) {
  edges <- results$edges[[c]] |> subset(!is_buff)
  
  data.frame(
    iso3c = c,
    country = iso3_codes[c],
    perc_inc = with(edges, sum(Ijk * cost_kmh) / sum(Ijk_orig * cost_kmh) - 1) * 100,
    cost = sum((edges$Ijk - edges$Ijk_orig) * edges$cost_kmh) / 1e6 * 10 # in billion of 2015 USD
  )
}) |> rbindlist() |> join(
  fread("results/grid_network/ALL_Mrealloc_fixed_irs_sigma38_rho0.csv") |> 
  compute(adj = TargetPerc / WgainPerc, keep = "iso3c"),
    on = "iso3c"
  ) |> 
  mutate(perc_inc_adj = perc_inc * adj,
         cost_adj = cost * adj) |> 
  with(set_names(cost_adj, country)) |> round(2)

# Heuristic calcualtions: same ratio, upgrading half as expensive
round(costs / tab$exp_match_realloc_perc * tab$ug_match_realloc_perc / 2, 2)


# Jobs gains
jobs_gains <- lapply(names(iso3_codes), function(c) {
  nodes <- results$nodes[[c]] |> subset(!is_buff)
  edges <- results$edges[[c]] |> subset(!is_buff)
  
  I_orig <- fmean(edges$Ijk_orig, w = nodes$pop_cell[edges$from] * nodes$pop_cell[edges$to])
  I_new <- fmean(edges$Ijk, w = nodes$pop_cell[edges$from] * nodes$pop_cell[edges$to]) 

  data.frame(
    country = iso3_codes[c],
    nodes = nrow(nodes),
    edges = nrow(edges),
    I_orig = I_orig,
    I_new = I_new,
    I_pch = (I_new / I_orig - 1) * 100,
    jobs_gains = -1/((1-0.33)*(1-I_orig)) * (I_new / I_orig - 1) * 100,
    markup_gains = -1/(1 - I_orig) * (I_new / I_orig - 1) * 100
  )
}) |> rbindlist()

jobs_gains |> 
  xtable::xtable(digits = 3) |> print(include.r = FALSE, booktabs = TRUE)
