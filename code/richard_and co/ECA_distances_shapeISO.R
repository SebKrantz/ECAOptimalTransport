library(fastverse)
fastverse_extend(sf, haven)
set_diag <- `diag<-`

ECA_centroids <- fread("data/ECA_shp_new/ECA_centroids.csv") |> rm_stub("_") |> 
  fsubset(is.finite(CX) & is.finite(CY)) |>
  collap( ~ countrycode + shapeISO)
# -> Matias fixed Kazhakstan, Armenia, and Georgia

ECA_centroids |> 
  fselect(shapeISO, shapeISO_nm, countrycode, lon = CX, lat = CY) |>
  fwrite("data/ECA_centroids.csv")

ECA_centroids %<>% st_as_sf(coords = c("CX", "CY"), crs = 4326)

# plot(ECA_centroids[, "shapeISO"])
# mapview::mapview(ECA_centroids)

# Now Computing Distance Matrix
library(osrm)
split_large_dist_matrix <- function(data, chunk_size = 100, verbose = FALSE) {
  n = nrow(data)
  res_list = list()
  
  # Loop over each chunk to compute the pairwise distances and travel times
  count = 0
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      # Define the row indices for the current chunks
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Extract the data for the current chunks
      ds_i = data[rows_i, ]
      ds_j = data[rows_j, ]
      
      if(verbose) {
        count = count + 1L
        cat(count," ")
      }
      
      # Perform the API call for the current chunks
      r_ij = osrmTable(src = ds_i, dst = ds_j, measure = c('duration', 'distance'))
      
      # Store the result in a list for later combination
      res_list[[paste(i, j, sep = "_")]] = r_ij
    }
  }
  
  # Combine the results from the list into one large matrix for durations and distances
  res_sources = matrix(NA, n, 2)
  res_destinations = matrix(NA, n, 2)
  res_durations = matrix(NA, n, n)
  res_distances = matrix(NA, n, n)
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Retrieve the result from the list
      r_ij = res_list[[paste(i, j, sep = "_")]]
      
      # Place the result into the corresponding location in the matrix
      res_sources[rows_i, ] = qM(r_ij$sources)
      res_destinations[rows_j, ] = qM(r_ij$destinations)
      res_durations[rows_i, rows_j] = r_ij$durations
      res_distances[rows_i, rows_j] = r_ij$distances
    }
  }
  
  # Create a result list to return
  res = list(
    sources = qDF(copyAttrib(mctl(res_sources), data)),
    destinations = qDF(copyAttrib(mctl(res_destinations), data)),
    durations = res_durations,
    distances = res_distances
  )
  
  rn = rownames(res$sources)
  if(length(rn) && suppressWarnings(!identical(as.integer(rn), seq_along(rn)))) {
    dimnames(res$durations) <- dimnames(res$distances) <- list(rn, rn)
  }
  
  return(res)
}

centroids <- ECA_centroids |> 
  ftransform(st_coordinates(geometry) |> mctl() |> set_names(c("lon", "lat"))) |>
  unclass() |> fselect(shapeISO, lon, lat) |> qM(1) |> qDF()
  
dist <- split_large_dist_matrix(centroids)
anyNA(dist$distances)
anyNA(dist$durations)
dist$destinations <- NULL
dist$centroids <- centroids
setrename(dist, sources = starts)

# Road Transport Cost Estimation
# ChatGPT suggests: Cost [€] = 1.05 × (distance in km) + 0.60 × (travel time in minutes)
# Iimi (2023) suggests: log(U.S. cents per ton-km) = 4.650 - 0.395 × log(speed in km/h) - 0.064 × log(distance in km) + 0.024 × crossborder dummy
# Also sector-specific results and different equations for domestic/international shipments...

# Compute travel cost
dist$cents_per_ton_km <- exp(4.650 - 0.395 * log(dist$distances / dist$durations * 60 / 1000) - 
                             0.064 * log(dist$distances / 1000) + 
                             0.024 * outer(ECA_centroids$countrycode, ECA_centroids$countrycode, "!="))
diag(dist$cents_per_ton_km) <- 0
dist$dollars_per_ton <- dist$cents_per_ton_km * dist$distances / 1e5

descr(vec(set_diag(dist$cents_per_ton_km, NA)))
descr(vec(set_diag(dist$dollars_per_ton, NA)))

# Counterfactual 1: Upgrading all roads to allow speeds of 70 km/h

descr(vec(dist$distances / dist$durations * 60 / 1000))
descr(vec(pmax(70, dist$distances / dist$durations * 60 / 1000)))

dist$cents_per_ton_km_g70kmh <- exp(4.650 - 0.395 * log(pmax(70, dist$distances / dist$durations * 60 / 1000)) - 
                                    0.064 * log(dist$distances / 1000) + 
                                    0.024 * outer(ECA_centroids$countrycode, ECA_centroids$countrycode, "!="))
diag(dist$cents_per_ton_km_g70kmh) <- 0
dist$dollars_per_ton_g70kmh <- dist$cents_per_ton_km_g70kmh * dist$distances / 1e5

descr(vec(set_diag(dist$cents_per_ton_km_g70kmh, NA)))
descr(vec(set_diag(dist$dollars_per_ton_g70kmh, NA)))

# Counterfactural 2: Eliminating border frictions

dist$cents_per_ton_km_no_border <- exp(4.650 - 0.395 * log(dist$distances / dist$durations * 60 / 1000) - 
                                       0.064 * log(dist$distances / 1000))
diag(dist$cents_per_ton_km_no_border) <- 0
dist$dollars_per_ton_no_border <- dist$cents_per_ton_km_no_border * dist$distances / 1e5

descr(vec(set_diag(dist$cents_per_ton_km_no_border, NA)))
descr(vec(set_diag(dist$dollars_per_ton_no_border, NA)))

saveRDS(dist, "data/ECA_centroids_distances_and_costs.rds")

# Combining results
result <- dist |>
  atomic_elem() |>
  unlist2d("variable", "from") |>
  pivot("from", names = list(from = "variable", to = "to"), how = "r") |>
  fsubset(from != to) |>
  fmutate(distances = distances / 1000) |>
  frename(shapeISO_o = from,
          shapeISO_d = to,
          distance_km = distances,
          duration_min = durations)

fnobs(result)

result |> fwrite("data/ECA_centroids_distances_and_costs.csv")

# Test joining to database

ECA <- read_dta("data/ECA_database_shapeISO.dta")
ECA |> join(result, on = c("shapeISO" = "shapeISO_o")) |> invisible() 
ECA |> join(result, on = c("shapeISO" = "shapeISO_d")) |> invisible() 



# Further counterfactuals: Adding new links -----------------------------------------

ECA_centroids <- fread("data/ECA_centroids.csv")
ECA_centroids %<>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
dist <- readRDS("data/ECA_centroids_distances_and_costs.rds")
net <- new.env()
load("data/ECA_centroids_network/ECA_centroids_network.RData", envir = net)
add_links <- net$add_links
edges <- net$edges
net <- sfnetworks::as_sfnetwork(edges, directed = FALSE)
nodes <- st_as_sf(net, "nodes")

# Check whether nodes align and computing correspondence for extended network
if(st_geometry(net, "nodes") |> st_distance(nodes, by_element = TRUE) |> allv(0) |> not()) stop("Mismatch")

# Match centroids to nodes
ECA_centroids$node <- dapply(st_distance(nodes, ECA_centroids), which.min)

# Check Ratios with Simulations
ratios <- sfnetworks::st_network_cost(net, weights = "duration", direction = "all")[ECA_centroids$node, ECA_centroids$node] / dist$durations
descr(vec(ratios))
ratios[!is.finite(ratios) | ratios < 1] <- 1
descr(vec(ratios))

# Now counterfactual: add links
edges |> with(distance / sp_distance) |> mean()
edges |> with(duration / sp_distance) |> mean()
edges |> with((distance / 1000) / (duration / 60)) |> descr()

add_links %<>% fmutate(sp_distance = st_length(geometry), 
                       distance = 1.444168 * sp_distance,
                       duration = (distance/1000) / 90 * 60) # 90kmh # Same as Existing: 0.00174973 * sp_distance

descr(add_links$duration)
add_links |> with((distance / 1000) / (duration / 60)) |> descr()
descr(edges$duration)
sum(add_links$distance) / sum(edges$distance)

# Creating extended network
net_ext <- sfnetworks::as_sfnetwork(rbind(fselect(edges, -passes), 
                                          fselect(add_links, -id)), directed = FALSE)

# Get index for matching
ind <- fmatch(mctl(st_coordinates(nodes)), mctl(st_coordinates(st_as_sf(net_ext, "nodes")))) # |> na_rm()

# New distances and durations
new_distances <- sfnetworks::st_network_cost(net_ext, weights = "distance", direction = "all")[ind, ind]
new_durations <- sfnetworks::st_network_cost(net_ext, weights = "duration", direction = "all")[ind, ind]

# Ratios
new_distance_ratios <- new_distances / sfnetworks::st_network_cost(net, weights = "distance", direction = "all")
descr(vec(new_distance_ratios))
new_distance_ratios[!is.finite(new_distance_ratios)] <- 1
new_distance_ratios[new_distance_ratios < 0.1] <- 0.1
descr(vec(new_distance_ratios))

new_duration_ratios <- new_durations / sfnetworks::st_network_cost(net, weights = "duration", direction = "all")
descr(vec(new_duration_ratios))
new_duration_ratios[!is.finite(new_duration_ratios)] <- 1
new_duration_ratios[new_duration_ratios < 0.1] <- 0.1
descr(vec(new_duration_ratios))

# Now Applying counterfactuals
dist$distances_net_ext <- dist$distances * new_distance_ratios[ECA_centroids$node, ECA_centroids$node]
dist$durations_net_ext <- dist$durations * new_duration_ratios[ECA_centroids$node, ECA_centroids$node]

# Compute travel cost
dist$cents_net_ext_per_ton_km <- exp(4.650 - 0.395 * log(dist$distances_net_ext / dist$durations_net_ext * 60 / 1000) - 
                                     0.064 * log(dist$distances_net_ext / 1000) + 
                                     0.024 * outer(ECA_centroids$countrycode, ECA_centroids$countrycode, "!="))
diag(dist$cents_net_ext_per_ton_km) <- 0
dist$dollars_net_ext_per_ton <- dist$cents_net_ext_per_ton_km * dist$distances_net_ext / 1e5

descr(vec(set_diag(dist$cents_net_ext_per_ton_km, NA)))
descr(vec(set_diag(dist$dollars_net_ext_per_ton, NA)))

# No border frictions

dist$cents_net_ext_per_ton_km_no_border <- exp(4.650 - 0.395 * log(dist$distances_net_ext / dist$durations_net_ext * 60 / 1000) - 
                                               0.064 * log(dist$distances_net_ext / 1000))
diag(dist$cents_net_ext_per_ton_km_no_border) <- 0
dist$dollars_net_ext_per_ton_no_border <- dist$cents_net_ext_per_ton_km_no_border * dist$distances_net_ext / 1e5

descr(vec(set_diag(dist$cents_net_ext_per_ton_km_no_border, NA)))
descr(vec(set_diag(dist$dollars_net_ext_per_ton_no_border, NA)))

saveRDS(dist, "data/ECA_centroids_distances_and_costs.rds")


# Middle corridor --------------------------------------------------------------------------

load("data/ECA_centroids_network/ECA_centroids_network.RData")
settfm(nodes, city_ctry = iif(is.na(city_ascii), paste0("N", seq_along(city)), paste(city_ascii, iso3, sep = " - ")))
if(any_duplicated(na_rm(nodes$city_ctry))) stop("Duplicates detected!")
edges_real <- qs::qread("data/ECA_centroids_network/edges_real_simplified.qs")
ECA_centroids <- fread("data/ECA_centroids.csv")
ECA_centroids %<>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
# Match centroids to nodes
ECA_centroids$node <- dapply(st_distance(nodes, ECA_centroids), which.min)

fastverse_extend(mapview, sfnetworks, tmap)

if(FALSE) mapview(nodes[, c("city_ctry", "population")], zcol = "city_ctry",
        legend = FALSE, map.types = c("Esri.WorldStreetMap", mapviewGetOption("basemaps"))) +
  mapview(add_links) +
  mapview(edges_real[, "passes"])

mc_nodes <- c("Tashkent - UZB", "Shymkent - KAZ", "Almaty - KAZ", "N828", # Khorgas
              "Usharal - KAZ", "N730", "Shalqar - KAZ", "N547->", "Beyneu - KAZ", "Aqtau - KAZ->", "Baku - AZE",
              "Shymkent - KAZ", "Shalqar - KAZ",
              "Tashkent - UZB", "Zhetisay - KAZ", "Beyneu - KAZ",
              "Baku - AZE", "Tbilisi - GEO", "Luleburgaz - TUR", # Istanbul
              "Luleburgaz - TUR", "Budapest - HUN"
)

# Add Node for Baku:
baku <- st_as_sfc(list(st_point(c(49.81750, 40.43022))), crs = 4326)
are <- mean(edges$sp_distance / edges$distance)
# Ferry connection from Aqtau - KAZ to Baku
aqtau_baku <- st_combine(c(nodes$geometry[nodes$city_ctry == "Aqtau - KAZ"], baku)) |> 
  st_cast("LINESTRING") |> st_as_sf() |> 
  fmutate(sp_distance = st_length(x), 
          distance = sp_distance,
          duration = distance / 1000 * 60/40) # 40 km/h
# Connect Baku to Tbilisi - GEO
baku_tbilisi <- st_combine(c(baku, nodes$geometry[nodes$city_ctry == "Tbilisi - GEO"])) |> 
  st_cast("LINESTRING") |> st_as_sf() |>  # |> mapview()
  fmutate(sp_distance = st_length(x), 
          distance = sp_distance / are, 
          duration = distance / 1000 * 60/90)
# Add direct connection from Shalqar (N547) to Beyneu
N547_beyneu <- st_combine(c(nodes$geometry[nodes$city_ctry == "N547"], 
                            nodes$geometry[nodes$city_ctry == "Beyneu - KAZ"])) |> 
  st_cast("LINESTRING") |> st_as_sf() |> # |> mapview()
  fmutate(sp_distance = st_length(x), 
          distance = sp_distance / are,
          duration = distance / 1000 * 60/90)
N604_N559 <- st_combine(c(nodes$geometry[nodes$city_ctry == "N604"], 
                          nodes$geometry[nodes$city_ctry == "N559"])) |> 
  st_cast("LINESTRING") |> st_as_sf() |> # |> mapview()
  fmutate(sp_distance = st_length(x), 
          distance = sp_distance / are,
          duration = distance / 1000 * 60/90)

# Test: with(N547_beyneu, (distance/1000) / (duration/60))

# Extended Network: 
edges_ext <- edges |> fselect(geometry, sp_distance, distance, duration) |> unclass() |> 
  Map(f = c, rowbind(aqtau_baku, baku_tbilisi, N547_beyneu, N604_N559)) |> qDF() |> st_as_sf()
net_ext <- as_sfnetwork(edges_ext, directed = FALSE)

# Check whether nodes align and computing correspondence for extended network
if(st_geometry(net, "nodes") |> st_distance(nodes, by_element = TRUE) |> allv(0) |> not()) stop("Mismatch")

# Get index for matching
ind <- fmatch(mctl(round(st_coordinates(nodes), 4)), mctl(round(st_coordinates(st_as_sf(net_ext, "nodes")), 4))) # |> na_rm()
fndistinct(ind)

# Now Computing MC Paths
tashkent_khorgas <- st_network_paths(net_ext, 
  from = ind[nodes$city_ctry == "Tashkent - UZB"], 
  to = ind[nodes$city_ctry == "N828"], weights = "duration", mode = "all")
usharal_shalkar <- st_network_paths(net_ext, 
  from = ind[nodes$city_ctry == "Usharal - KAZ"], 
  to = ind[nodes$city_ctry == "Shalqar - KAZ"], weights = "duration", mode = "all")
tashkent_shalkar <- st_network_paths(net_ext, 
  from = ind[nodes$city_ctry == "Tashkent - UZB"], 
  to = ind[nodes$city_ctry == "Shalqar - KAZ"], weights = "duration", mode = "all")
tashkent_beyneu <- st_network_paths(net_ext, # Feeder route
  from = ind[nodes$city_ctry == "Tashkent - UZB"], 
  to = ind[nodes$city_ctry == "Beyneu - KAZ"], weights = "duration", mode = "all")
shalkar_baku <- st_network_paths(net_ext, 
  from = ind[nodes$city_ctry == "Shalqar - KAZ"], 
  to = which.min(st_distance(st_as_sf(net_ext, "nodes"), baku)), weights = "duration", mode = "all")
baku_istanbul <- st_network_paths(net_ext, 
  from = which.min(st_distance(st_as_sf(net_ext, "nodes"), baku)), 
  to = ind[nodes$city_ctry == "Luleburgaz - TUR"], weights = "duration", mode = "all")
istanbul_budapest <- st_network_paths(net_ext, 
  from = ind[nodes$city_ctry == "Luleburgaz - TUR"], 
  to = ind[nodes$city_ctry == "Budapest - HUN"], weights = "duration", mode = "all")

mc_paths <- rowbind(tashkent_khorgas, tashkent_shalkar, tashkent_beyneu, usharal_shalkar, shalkar_baku, 
                    baku_istanbul, istanbul_budapest)
mc_paths <- list(nodes = unique(unlist(mc_paths$node_paths)), 
                 edges = unique(unlist(mc_paths$edge_paths)))

# Ensuring speed is 90km/h along corridor routes
edges_ext <- rowbind(edges_ext |> ss(1:(nrow(edges_ext)-4)) |> 
                       fmutate(duration = iif(seq_along(duration) %in% mc_paths$edges, distance / 1000 / 90 * 60, duration)), 
                     edges_ext |> ss((nrow(edges_ext)-3):nrow(edges_ext)))

net_ext <- as_sfnetwork(edges_ext, directed = FALSE)
nodes_ext <- st_geometry(net_ext, "nodes") |> st_as_sf()
edges_real_ext <- edges_real |> fselect(geometry, sp_distance, distance, duration, passes) |> unclass() |> 
  Map(f = c, rowbind(aqtau_baku, baku_tbilisi, N547_beyneu, N604_N559) |> fmutate(passes = 1)) |> qDF() |> st_as_sf()
tfm(edges_real_ext) <- qDF(edges_ext) |> fselect(distance, duration)
settfm(edges_real_ext, speed_kmh = distance / duration * 60 / 1000)

tmap_options(raster.max_cells = 1e8)
tmap_mode("plot")

pdf("figures/ECA_centroids_network_actual_discretized_middle_corridor.pdf", width = 17, height = 7)
tm_basemap("CartoDB.Positron", zoom = 5) +
  tm_shape(ss(edges_real_ext, -mc_paths$edges)) +
  tm_lines(col = "grey70", lwd = 2) +
  tm_shape(ss(edges_real_ext, mc_paths$edges)) +
  tm_lines(col = "orange", lwd = 2) +
  tm_shape(ss(nodes_ext, -mc_paths$nodes)) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(ss(nodes_ext, mc_paths$nodes)) + tm_dots(size = 0.2) +
  # tm_shape(subset(nodes, mc_nodes)) + tm_dots(size = 0.2, fill = "red") +
  tm_layout(frame = FALSE) 
dev.off()

pdf("figures/ECA_centroids_network_actual_discretized_middle_corridor_speed_EWTM.pdf", width = 17, height = 7)
tm_basemap("Esri.WorldTopoMap", zoom = 6) +
  tm_shape(edges_real_ext) +
  tm_lines(col = "speed_kmh",
           col.scale = tm_scale_continuous(5, values = "turbo"),
           col.legend = tm_legend("Speed in km/h", position = c("left", "top"), frame = FALSE, 
                                  text.size = 1, title.size = 1.2, margins = c(0, -0.5, 0, 0),
                                  title.padding = c(0, 0, -0.5, 0), 
                                  item.space = 0, item.height = 2, item.width = 0.5)) +
  tm_shape(subset(nodes, key_city)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !key_city)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))
dev.off()

# New distances and durations
new_distances <- sfnetworks::st_network_cost(net_ext, weights = "distance", direction = "all")[ind, ind]
new_durations <- sfnetworks::st_network_cost(net_ext, weights = "duration", direction = "all")[ind, ind]

# Ratios
new_distance_ratios <- new_distances / sfnetworks::st_network_cost(net, weights = "distance", direction = "all")
descr(vec(new_distance_ratios))
new_distance_ratios[!is.finite(new_distance_ratios)] <- 1
new_distance_ratios[new_distance_ratios < 0.1] <- 0.1
descr(vec(new_distance_ratios))

new_duration_ratios <- new_durations / sfnetworks::st_network_cost(net, weights = "duration", direction = "all")
descr(vec(new_duration_ratios))
new_duration_ratios[!is.finite(new_duration_ratios) | new_duration_ratios > 1] <- 1
new_duration_ratios[new_duration_ratios < 0.1] <- 0.1
descr(vec(new_duration_ratios))

# Now Applying counterfactuals
dist$distances_mc <- dist$distances * new_distance_ratios[ECA_centroids$node, ECA_centroids$node]
dist$durations_mc <- dist$durations * new_duration_ratios[ECA_centroids$node, ECA_centroids$node]

# Compute travel cost
dist$cents_mc_per_ton_km <- exp(4.650 - 0.395 * log(dist$distances_mc / dist$durations_mc * 60 / 1000) - 
                                0.064 * log(dist$distances_mc / 1000) + 
                                0.024 * outer(ECA_centroids$countrycode, ECA_centroids$countrycode, "!="))
diag(dist$cents_mc_per_ton_km) <- 0
dist$dollars_mc_per_ton <- dist$cents_mc_per_ton_km * dist$distances_mc / 1e5

descr(vec(set_diag(dist$cents_mc_per_ton_km, NA)))
descr(vec(set_diag(dist$dollars_mc_per_ton, NA)))

# No border frictions

dist$cents_mc_per_ton_km_no_border <- exp(4.650 - 0.395 * log(dist$distances_mc / dist$durations_mc * 60 / 1000) - 
                                          0.064 * log(dist$distances_mc / 1000))
diag(dist$cents_mc_per_ton_km_no_border) <- 0
dist$dollars_mc_per_ton_no_border <- dist$cents_mc_per_ton_km_no_border * dist$distances_mc / 1e5

descr(vec(set_diag(dist$cents_mc_per_ton_km_no_border, NA)))
descr(vec(set_diag(dist$dollars_mc_per_ton_no_border, NA)))

saveRDS(dist, "data/ECA_centroids_distances_and_costs.rds")


# Combining results ----------------------------


result <- dist |>
  atomic_elem() |>
  unlist2d("variable", "from") |>
  pivot("from", names = list(from = "variable", to = "to"), how = "r") |>
  fsubset(from != to) |>
  ftransformv(c(distances, distances_net_ext, distances_mc), `/`, 1000) |>
  frename(shapeISO_o = from,
          shapeISO_d = to,
          distance_km = distances, distance_net_ext_km = distances_net_ext, distance_mc_km = distances_mc,
          duration_min = durations, duration_net_ext_min = durations_net_ext, duration_mc_min = durations_mc)

fnobs(result)
qsu(result)

result |> fwrite("data/ECA_centroids_distances_and_costs.csv")

