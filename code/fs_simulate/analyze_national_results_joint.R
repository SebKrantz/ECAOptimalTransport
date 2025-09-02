#####################################################################
# Transport Network: Analyze Optimal Trans-African GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, units, sfnetworks, tmap, ggplot2, install = TRUE)
fastverse_conflicts()

files <- list.files("results/grid_network/country", pattern = "_sigma2_") |> 
  grep(pattern = "_irs_", invert = TRUE, value = TRUE) |> 
  grep(pattern = "_ann_", value = TRUE, invert = TRUE) |>
  grep(pattern = "_ug_", value = TRUE, invert = TRUE)
countries <- substr(files, 15, 17) |> unique()
setdiff(unique(substr(list.files("data/grid_network/country"), 1, 3)), countries)

results <- list(
  nodes = lapply(files[startsWith(files, "nodes")], function(f) fread(paste0("results/grid_network/country/", f))) |>
          set_names(substr(files[startsWith(files, "nodes")], 15, 17)),
  edges = lapply(files[startsWith(files, "edges")], function(f) fread(paste0("results/grid_network/country/", f))) |>
          set_names(substr(files[startsWith(files, "edges")], 15, 17))
)

# results %<>% rapply2d(function(x) subset(x, is_buff))

# Statistics on the upgrade extent -----------------------------------------------------------------

# Cost of all possible work (should be 13.5 billion in millions)
(cost <- results$edges |> sapply(with, sum(cost_kmh * time_efficiency)))
# Budget spent (should give 1 or 2 billion in millions)
(bud <- results$edges |> sapply(with, sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 70)-Ijk_orig), 0), 0)*cost_kmh*time_efficiency)))
# Amount spent on different types of work (1 = new construction, 0 = upgrade)
sort(bud / cost * 100)
# TODO: should be 10%

# Road km built/upgraded
results$edges |> sapply(with, sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 70)-Ijk_orig), 0), 0)*sp_distance))

# Number of roads worked on (TRUE) (extensive margin)
results$edges |> sapply(with, table(Ijk-Ijk_orig > 1)) # |> proportions()

# Work Intensity (intensive margin = km/h added)
results$edges |> sapply(with, fmedian((Ijk-Ijk_orig)[Ijk-Ijk_orig > 1])) |> sort()

# Examining percentage on infra links on edges between urban centers -----------------------------------------------------------------

# Ana: When discussing the reasons behind sub-optimal investments, I would like us to think if we can come with a systematic measure 
# of the percentage of peripheral investments relative to intercity ones. Visually the root of the problem is clear. But we need to think 
# about formalizing it more. I have a couple of ideas. Let’s talk next week. 

compute_cent_edges <- function(nodes, edges) {
  
  nodes$big_city <- nodes$product > 5 # $pop_cell > quantile(nodes$pop_cell, 0.75)
  
  g <- cppRouting::makegraph(edges |> select(from, to, cost = duration))
  mp <- expand.grid(o = which(nodes$big_city), d = which(nodes$big_city)) |> subset(o != d)
  paths <- cppRouting::get_multi_paths(g, mp$o, mp$d, long = TRUE) |> unique()
  paths$id <- group(paths$from, paths$to)
  paths <- roworder(paths, id)
  paths$from <- paths$node
  paths$to <- flag(paths$node, -1, g = paths$id)
  paths <- na_omit(paths)
  
  edg <- unique(fmatch(select(paths, from, to), select(edges, from, to)))
  edg2 <- unique(fmatch(select(paths, to, from), select(edges, from, to)))
  unique(na_rm(c(edg, edg2))) |> as.integer() |> sort()
}

compute_cent_edges_perc <- function(nodes, edges) {
  
  mpe <- compute_cent_edges(nodes, edges)
  
  tot_realloc <- (edges[["Ijk"]] - edges[["Ijk_orig"]]) * edges$distance
  
  # This is the percentage devoted to important links.
  200 * sum(abs(tot_realloc[mpe])) / sum(abs(tot_realloc))
  
}

cent_edges_perc <- sapply(countries, function(c) {

  nodes <- results$nodes[[c]]
  edges <- results$edges[[c]]
  
  compute_cent_edges_perc(nodes, edges)

  
})

cent_edges_perc |> round(2) |> sort()

# cent_edges_results <- list()
cent_edges_results$`10perc_fixed` <- cent_edges_perc
cent_edges_results_df <- unlist2d(cent_edges_results, "spec") |> 
  transpose(make.names = 1, keep.names = "country")

writexl::write_xlsx(cent_edges_results_df, "results/grid_network/ALL_cent_edges_perc.xlsx")


# Statistics on the economic gains -----------------------------------------------------------------

# Global Welfare Gains (Ratio)
wgains <- results$nodes |> lapply(subset, !is_buff) |> 
  sapply(with, sum(uj * Lj_orig) / sum(uj_orig * Lj_orig)) |> 
  subtract(1) |> multiply_by(100) |> round(5)

# Test
wgains2 <- results$nodes |> lapply(subset, !is_buff) |> 
  sapply(with, sum((Cj/Lj/0.4)^0.4 * Lj_orig) / sum((Cj_orig/Lj_orig/0.4)^0.4 * Lj_orig)) |> 
  subtract(1) |> multiply_by(100) |> round(5)

all.equal(wgains, wgains2)

# Consumption Gains (Ratio)
cgains <- results$nodes |> lapply(subset, !is_buff) |> 
  sapply(with, sum(Cj) / sum(Cj_orig))

# Now applying consumption gains proportionally to all locations
welfare <- results$nodes |> lapply(subset, !is_buff) |> lapply(with,
    list(welfare = sum(uj * Lj_orig),
         welfare_orig = sum(uj_orig * Lj_orig),
         welfare_cf = sum(((sum(Cj)/sum(Cj_orig)*Cj_orig)/Lj/0.4)^0.4 * Lj_orig))                                                         
)

# Computing decomposition 
welfare <- welfare |> lapply(with, 
  list("Total welfare change (%)" = (welfare - welfare_orig) / welfare_orig * 100,
       "Volume effect (%)" = (welfare_cf - welfare_orig) / welfare_orig * 100,
       "Allocation effect (%)" = (welfare - welfare_cf) / welfare_orig * 100)
) |> rowbind(idcol = "iso3c", return = "data.frame")

# Test
all.equal(unattrib(((cgains)^0.4-1) * 100), welfare$`Volume effect (%)`)

# Correlates
if(!exists("predictors")) {
  predictors <- africamonitor::am_data(ctry = names(wgains), series = c("NY_GDP_PCAP_KD", "AG_LND_TOTL_K2", "LP_LPI_OVRL_XQ", "LP_LPI_INFR_XQ", "IC_BUS_EASE_DFRN_XQ_DB1719")) |> 
    collap(~ISO3, flast, flast) 
}
predictors %<>% subset(match(names(wgains), ISO3))
rbind(pearson = drop(pwcor(wgains, nv(predictors))),
      spearman = drop(pwcor(wgains, nv(predictors), method = "spearman")),
      kendall = drop(pwcor(wgains, nv(predictors), method = "kendall")))


barplot(sort(wgains), horiz = TRUE, las = 2)

wgains_df <- janitor::clean_names(welfare) |> 
  rename(total_welfare_change_percent = wgain) |> 
  roworder(wgain) |> 
  mutate(country = countrycode::countrycode(iso3c, "iso3c", "country.name.en") |> 
           replace_na("Kosovo") |> qF(sort = FALSE))

# Export 
wgains_df |> 
  colorder(iso3c, country) |> 
  writexl::write_xlsx("results/grid_network/welfare_gains_sigma2_irs.xlsx")

# Plot
wgains_df |> 
  ggplot(aes(x = wgain, y = country, fill = wgain)) +
    geom_bar(stat = "identity") +
    geom_bar(aes(x = volume_effect_percent), fill = "white", alpha = 0.3, stat = "identity") +
    geom_text(aes(label = paste0(signif(wgain, 2), "%")), nudge_x = 0.02, size = 3) +
    viridis::scale_fill_viridis() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                       expand = expansion(add = c(0, 0.03))) +
    theme_minimal() + guides(fill = "none") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.ticks.length.x = unit(0.4, "lines"),
          axis.ticks.length.y = unit(0, "lines")) +
    labs(x = "Welfare Gain from Optimal Network Reallocation", y = "Country")
  
ggsave("figures/ECA_grid_network_realloc_barchart.pdf", width = 7, height = 7)

# Joint Plot
wgains_df %<>% join(wgains_df_drs, on = "country", suffix = "_drs", how = "inner")

# Reallocation Plot
wgains_df |> 
  ggplot(aes(x = wgain, y = country)) + # , fill = wgain
  geom_bar(stat = "identity", fill = "orange") +
  geom_bar(aes(x = volume_effect_percent), fill = "white", alpha = 0.3, stat = "identity") +
  geom_point(aes(x = wgain_drs), shape = 108, size = 5) +
  geom_point(aes(x = volume_effect_percent_drs), shape = 108, size = 5, colour = "grey50") +
  geom_text(aes(x = wgain_drs, label = paste0(round(wgain_drs, 2), "%")), nudge_x = 0.028, size = 3) +
  geom_text(aes(label = paste0(round(wgain, 2), "%")), nudge_x = 0.028 + abs(0.17-wgains_df$wgain)/4*(wgains_df$wgain < 0.17), size = 3) +
  annotate(
    geom = "text", x = 0.5, y = 10.5, 
    label = expression(bold("\u2014")~" Decreasing Returns: "~gamma < beta), size = 4
  ) +
  annotate(
    geom = "text", x = 0.5, y = 6, 
    label = expression("Lighter Color: Volume Effect\n(Total = Volume + Allocation)"), size = 4
  ) +
  # viridis::scale_fill_viridis(option = "F", begin = 0.4, end = 1, direction = -1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 0.03)), limits = c(0, 0.72)) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Welfare Gain from Optimal Network Reallocation", y = "Country")

ggsave("figures/ECA_grid_network_realloc_barchart_both_decomp_nograd.pdf", width = 7, height = 7)

# 10% Expansion Plot
wgains_df |> 
  ggplot(aes(x = wgain, y = country, fill = wgain)) +
  geom_bar(stat = "identity") +
  geom_point(aes(x = wgain_drs), shape = 108, size = 5) +
  geom_text(aes(x = wgain_drs, label = paste0(round(wgain_drs, 2), "%")), nudge_x = 0.02, size = 3) +
  geom_text(aes(label = paste0(round(wgain, 2), "%")), nudge_x = 0.02 + abs(0.1-wgains_df$wgain)/4*(wgains_df$wgain < 0.07), size = 3) +
  annotate(
    geom = "text", x = 0.3, y = 10, 
    label = expression(bold("\u2014")~" Decreasing Returns: "~gamma < beta), size = 4
  ) +
  viridis::scale_fill_viridis(option = "D", begin = 0.2, end = 1, direction = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 0.03)), limits = c(0, 0.38)) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Welfare Gain from Optimal 10% Network Expansion", y = "Country")

ggsave("figures/ECA_grid_network_10perc_barchart_both.pdf", width = 7, height = 7)


# 10% Expansion Plot: joint
wgains_df |> 
  ggplot(aes(x = wgain, y = country)) +
  geom_bar(stat = "identity", fill = "orange") +
  geom_bar(aes(x = volume_effect_percent), fill = "white", alpha = 0.3, stat = "identity") +
  geom_point(aes(x = wgain_drs), shape = 108, size = 5) +
  geom_point(aes(x = volume_effect_percent_drs), shape = 108, size = 5, colour = "grey50") +
  geom_text(aes(x = wgain_drs, label = paste0(round(wgain_drs, 2), "%")), nudge_x = 0.02, size = 3) +
  geom_text(aes(label = paste0(round(wgain, 2), "%")), nudge_x = 0.02 + abs(0.1-wgains_df$wgain)/4*(wgains_df$wgain < 0.07), size = 3) +
  annotate(
    geom = "text", x = 0.22, y = 10.5, 
    label = expression(bold("\u2014")~" Decreasing Returns: "~gamma < beta), size = 4
  ) +
  annotate(
    geom = "text", x = 0.22, y = 6, 
    label = expression("Lighter Color: Volume Effect\n(Total = Volume + Allocation)"), size = 4
  ) +
  # viridis::scale_fill_viridis(option = "D", begin = 0.3, end = 0.95, direction = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 0.03)), limits = c(0, 0.38)) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Welfare Gain from Optimal 10% Network Expansion", y = "Country")

ggsave("figures/ECA_grid_network_10perc_barchart_both_decomp_nograd.pdf", width = 7, height = 7)



# 10% Upgrading Plot
wgains_df |> 
  ggplot(aes(x = wgain, y = country)) + # , fill = wgain
  geom_bar(stat = "identity", fill = "orange") +
  geom_bar(aes(x = volume_effect_percent), fill = "white", alpha = 0.3, stat = "identity") +
  geom_point(aes(x = wgain_drs), shape = 108, size = 5) +
  geom_point(aes(x = volume_effect_percent_drs), shape = 108, size = 5, colour = "grey50") +
  geom_text(aes(x = wgain_drs, label = paste0(round(wgain_drs, 2), "%")), nudge_x = 0.014, size = 3) +
  geom_text(aes(label = paste0(round(wgain, 2), "%")), nudge_x = 0.014 + abs(0.07-wgains_df$wgain)/3*(wgains_df$wgain < 0.07), size = 3) +
  annotate(
    geom = "text", x = 0.12, y = 9, hjust = 0, vjust = 0.5,
    label = expression("Decreasing Returns: "~gamma < beta), size = 4
  ) +
  annotate("rect", xmin = 0.105, xmax = 0.114, ymin = 6.6, ymax = 7.6, fill = "orange", alpha = 0.7) +
  annotate("point", x = 0.11, y = 9.1, shape = 108, size = 6, colour = "black") +
  annotate("point", x = 0.11, y = 7.1, shape = 108, size = 6, colour = "grey50") +
  annotate(
    geom = "text", x = 0.12, y = 7, hjust = 0, vjust = 0.5,
    label = expression("Volume Effect (Total = Volume + Allocation)"), size = 4
  ) +
  # viridis::scale_fill_viridis(option = "H", begin = 0.3, end = 0.95, direction = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 0.03)), limits = c(0, 0.275)) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Welfare Gain from Optimal 10% Network Upgrading", y = "Country")

ggsave("figures/ECA_grid_network_upgrading_barchart_both_decomp_nograd_legadj.pdf", width = 7, height = 7)


# Excusus: Upgrading Amount to Match Rellocation Gains
Mrealloc <- list(DRS = fread("results/grid_network/ALL_Mrealloc_fixed_sigma38_rho0.csv"),
                 IRS = fread("results/grid_network/ALL_Mrealloc_fixed_irs_sigma38_rho0.csv"))
Mrealloc %<>% lapply(. %>% mutate(PercIncAdj = PercInc * TargetPerc / WgainPerc) %>%
                     roworder(PercIncAdj) %>% 
                     mutate(country = countrycode::countrycode(iso3c, "iso3c", "country.name.en") |> 
                                replace_na("Kosovo") |> qF(sort = FALSE)))
# Only DRS
Mrealloc$DRS |> 
  ggplot(aes(x = PercIncAdj, y = country)) + # , fill = wgain
  geom_bar(stat = "identity", fill = "orange") +
  geom_text(aes(label = paste0(round(PercIncAdj, 2), "%")), nudge_x = 2.5, size = 3) +
  # viridis::scale_fill_viridis(option = "F", begin = 0.4, end = 1, direction = -1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 7))) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Increase in Infrastructure Stock to Match Optimal Reallocation Gains", y = "Country")

ggsave("figures/ECA_grid_network_Mrealloc_barchart_DRS_nograd.pdf", width = 7, height = 7)



# Both 
Mrealloc$IRS |> 
  select(country, PercIncAdj) |> 
  join(Mrealloc$DRS |> select(country, PercIncAdj), on = "country", suffix = "_drs") |> 

  ggplot(aes(x = PercIncAdj, y = country)) + # , fill = wgain
  geom_bar(stat = "identity", fill = "orange") +
  geom_point(aes(x = PercIncAdj_drs), shape = 108, size = 5) +
  geom_text(aes(x = PercIncAdj_drs, label = paste0(round(PercIncAdj_drs, 2), "%")), nudge_x = 8, size = 3) +
  geom_text(aes(label = paste0(round(PercIncAdj, 2), "%")), nudge_x = 8, size = 3) +
  annotate(
    geom = "text", x = 100, y = 10.5, 
    label = expression(bold("\u2014")~" Decreasing Returns: "~gamma < beta), size = 4
  ) +
  # viridis::scale_fill_viridis(option = "F", begin = 0.4, end = 1, direction = -1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1, scale = 1), 
                     expand = expansion(add = c(0, 10))) +
  theme_minimal() + guides(fill = "none") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length.x = unit(0.4, "lines"),
        axis.ticks.length.y = unit(0, "lines")) +
  labs(x = "Increase in Infrastructure Stock to Match Optimal Reallocation Gains", y = "Country")

ggsave("figures/ECA_grid_network_Mrealloc_barchart_both_nograd.pdf", width = 7, height = 7)


# Tables for Export
eca_regions <- list(
  "Central Asia" = c("Kazakhstan", "Kyrgyz Republic", "Tajikistan", "Turkmenistan", "Uzbekistan"),
  "Central Europe" = c("Bulgaria", "Croatia", "Poland", "Romania"),
  "Eastern Europe" = c("Belarus", "Moldova", "Ukraine"),
  "Russian Federation" = c("Russian Federation"),
  "South Caucasus" = c("Armenia", "Azerbaijan", "Georgia"),
  "Türkiye" = c("Türkiye"),
  "Western Balkans" = c("Albania", "Bosnia and Herzegovina", "Kosovo", "Montenegro", "Republic of North Macedonia", "Serbia")
) %>% lapply(list) %>% rowbind(idcol = "cluster", use.names = FALSE) %>% 
  mutate(iso3c = countrycode::countryname(V1, "iso3c") |> replace_na("XKX"), 
         V1 = NULL)


result <- Mrealloc |> 
  lapply(function(x) select(x, country, iso3c, 
         realloc_gain_perc = TargetPerc, 
         exp_match_realloc_perc = PercIncAdj) |> 
         join(eca_regions, on = "iso3c") |> 
         colorder(cluster))

result |> 
  writexl::write_xlsx("results/grid_network/ALL_Mrealloc_table.xlsx")

# Excursus: Database of control characteristics at the country level

countries <- unique(substr(list.files("data/grid_network/country"), 1, 3))
data <- sapply(countries, function(x) {
  list(nodes = fread(sprintf("data/grid_network/country/%s_nodes.csv", x)),
       edges = fread(sprintf("data/grid_network/country/%s_edges.csv", x)))       
}, simplify = FALSE)

results <- data |> lapply(function(x) {
  x$nodes |> subset(!is_buff) |> 
    summarise(GDP = fsum(predicted_GCP_const_2017_PPP),
              POP = fsum(pop_cell),
              GDPCap = fmean(cell_GDPC_const_2017_PPP, w = pop_cell), 
              rugg = fmean(rugg), 
              rugg_wavg = fmean(rugg, w = pop_cell),
              pop_wpop_km2 = fmean(pop_wpop_km2), 
              cost_kmh = fmean(cost_kmh),
              cost_kmh_wavg = fmean(cost_kmh, w = pop_cell),
              own_product = fsum(own_product)) |> 
add_vars(
  x$edges |> mutate(gravity = x$nodes$pop_cell[from] * x$nodes$pop_cell[to] / sp_distance) |> 
    subset(!is_buff) |> 
    summarise(across(c(distance, duration, speed_kmh, route_efficiency, time_efficiency), fmean),
              across(c(distance, duration, speed_kmh, route_efficiency, time_efficiency), 
                     list(wavg = fmean), w = gravity, .names = TRUE))
)
}) |> rowbind(idcol = "iso3c")

setrelabel(results, 
  GDP = "GDP (PPP, 2017 constant US$)",
  POP = "Population (millions)",
  GDPCap = "GDP per Capita (PPP, 2017 constant US$)",
  rugg = "Ruggedness (mean)",
  rugg_wavg = "Ruggedness (weighted mean)",
  pop_wpop_km2 = "Population Density (people per km²)",
  cost_kmh = "Cost per km/h (mean)",
  cost_kmh_wavg = "Cost per km/h (weighted mean)",
  own_product = "Own Product (sum) (number of cities producing a unique product)",
  distance = "Distance (mean)",
  duration = "Duration (mean)",
  speed_kmh = "Speed (mean)",
  route_efficiency = "Route Efficiency (mean)",
  time_efficiency = "Time Efficiency (mean)",
  distance_wavg = "Distance (gravity weighted mean)",
  duration_wavg = "Duration (gravity weighted mean)",
  speed_kmh_wavg = "Speed (gravity weighted mean)",
  route_efficiency_wavg = "Route Efficiency (gravity weighted mean)",
  time_efficiency_wavg = "Time Efficiency (gravity weighted mean)"
)

# Add other data
results <- fread("results/grid_network/control_dataset.csv")

data <- am_data(ctry = results$iso3c, 
                series = c("AG_LND_TOTL_K2", "LP_LPI_OVRL_XQ", "LP_LPI_INFR_XQ",
                           "IC_BUS_EASE_DFRN_XQ_DB1719", "TRD_ACRS_BRDR_DB1619_DFRN")) |> 
        collap( ~ ISO3, flast, cols = is.numeric) 

results %<>% join(data, on = c("iso3c" = "ISO3"))

results |> fwrite("results/grid_network/control_dataset.csv")
results |> haven::write_dta("results/grid_network/control_dataset.dta")

# Adding consolidate data
cons_data <- haven::read_dta("misc/Consolidated and control data.dta")
cons_data |> join(results, on = c("iso3" = "iso3c"), drop = "x") |> View()
  haven::write_dta("results/grid_network/consolidated_and_control_dataset.dta")

# Country Plots --------------------------------------------------------------------------------------------

all_identical(lapply(results, names))
t_res <- t_list(results)
extra_countries = ""

save_plot = TRUE
outfolder = "figures/grid_network/10perc_fixed_ce"
qualifier = "10pf"

country_names <- countrycode::countrycode(names(t_res), "iso3c", "country.name.en") |> replace_na("Kosovo")
names(country_names) <- names(t_res)

for (country in names(t_res)) {
  
  .c(nodes, edges) %=% t_res[[country]]
  
  # Import
  # case_poly = polys[polys$NAME == country,]
  I_stat <- edges$Ijk_orig # read.csv(paste0("./Build/temp/speed/speed_", country, ".csv"), header = T)
  I_opt <- edges$Ijk # read.csv(paste0(infolder, "/Optimised_Networks/", country, ".csv"), header = F)
  I_diff <- I_opt - I_stat
  outcomes <- nodes # read.csv(paste0(infolder, "/Network_outcomes/", country, "_outcomes.csv"))
  
  # Scale
  # I_opt_scaled = 13*(sqrt(I_opt-4) / max(sqrt(I_opt)))
  # I_stat_scaled = 13*(sqrt(I_stat-4) / max(sqrt(I_opt)))
  
  if (nrow(nodes) > 350) {
    I_opt_scaled <- 5 * (I_opt^0.8 / max(I_opt^0.8))
    I_stat_scaled <- 5 * (I_stat^0.8 / max(I_opt^0.8))
    I_diff_scaled <- 5 * (I_diff^0.8 / max(I_diff^0.8))
  } else {
    I_opt_scaled <- 10 * (I_opt^0.8 / max(I_opt^0.8))
    I_stat_scaled <- 10 * (I_stat^0.8 / max(I_opt^0.8))
    I_diff_scaled <- 10 * (I_diff^0.8 / max(I_diff^0.8))
  }
  
  outcomes$pop_scaled <- 4 * (outcomes$pop_cell^0.33 / max(outcomes[!outcomes$is_buff, "pop_cell"])^0.33)
  outcomes[outcomes$is_buff, "color"] <- scales::alpha("grey", .5)
  
  outcomes$different_good <- outcomes$own_product
  outcomes$zeta <- outcomes$uj / outcomes$uj_orig
  outcomes$zetacol <- cut(outcomes$zeta, quantile(outcomes$zeta, seq(0, 1, 0.05))) # TG does it globally
  I_diff_q <- quantile(I_diff[!edges$is_buff], seq(0, 1, 1/10)) # seq(0, 1, 1/9)
  I_diff_col <- cut(I_diff[!edges$is_buff], I_diff_q, include.lowest = TRUE)
  
  for (graph in "cent") { # c("stat", "opt")) 
    # Scale
    if (graph == "opt") {
      edges[!edges$is_buff, "color"] <- cols4all::c4a("matplotlib.seismic", 16)[unclass(I_diff_col)+3L+sum(I_diff_q > 0)-5L] # RColorBrewer::brewer.pal(9, "YlOrRd")[I_diff_col]
      outcomes[!outcomes$is_buff, "color"] <- viridis::turbo(26, direction = 1)[unclass(outcomes$zetacol[!outcomes$is_buff])+3L]
    } else if (graph == "stat") {
      outcomes[!outcomes$is_buff, "color"] <- "dodgerblue4"
    } else { 
      outcomes[!outcomes$is_buff, "color"] <- "dodgerblue4"
      mpe <- compute_cent_edges(nodes, edges)
      edges[!edges$is_buff, "color"] <- "grey30"
      edges[mpe, "color"] <- "red3"
      I_cent_scaled <- I_opt_scaled
    }
    
    if(save_plot) pdf(file = paste(outfolder, "/", country, "_", graph, "_", qualifier, ".pdf", sep = ""), width = 15, height = 11)
    
    plot(outcomes$pwx, outcomes$pwy, main = country_names[country], bty = "n", 
         pch = ifelse(!outcomes$is_buff, 19, 1), axes = F, 
         ylab = "", xlab = "", asp = 1, type = "n") # if you want to plot it
    
    points(outcomes$pwx[outcomes$is_buff], outcomes$pwy[outcomes$is_buff], pch = 21, 
           bg = outcomes$color[outcomes$is_buff], 
           col = ifelse(outcomes$different_good[outcomes$is_buff] == 1, "black", NA), 
           cex = outcomes$pop_scaled[outcomes$is_buff])
    
    # plot(case_poly, add = T, col = alpha("lightskyblue1", 1), border = NA)
    
    for (i in seq_row(edges)) {
      points(c(edges$from_lon[i], edges$to_lon[i]), c(edges$from_lat[i], edges$to_lat[i]),
             type = "l", lwd = get(paste0("I_", graph, "_scaled"))[i], # I_diff_scaled[i]
             col = ifelse(edges$is_buff[i], scales::alpha("grey", .2), scales::alpha(edges$color[i], .8))) # scales::alpha("dodgerblue", .8)
    }
    
    # points(outcomes$pwx, outcomes$pwy, pch = ifelse(!outcomes$is_buff, 19, 19), col=outcomes$color, cex = outcomes$pop_scaled)
    points(outcomes$pwx[!outcomes$is_buff], outcomes$pwy[!outcomes$is_buff], 
           pch = 21, lwd = 2, 
           bg = outcomes$color[!outcomes$is_buff], 
           col = ifelse(outcomes$different_good[!outcomes$is_buff] == 1, "white", NA), 
           cex = outcomes$pop_scaled[!outcomes$is_buff])
    
    
    if(save_plot) dev.off()
  }
  
  print(country)
}

  
  
  

