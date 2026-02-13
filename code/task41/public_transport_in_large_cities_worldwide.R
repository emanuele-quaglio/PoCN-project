#!/usr/bin/env Rscript
# ==============================================================================
# public_transport_in_large_cities_worldwide.R
#
# Standalone pipeline merging the full workflow from:
#   1. inspection_reconstruction.rmd  (data inspection, JSON->Parquet, graph CSV)
#   2. analysis.rmd                   (network analytics, plotting, global plots)
#
# Usage:
#   Rscript public_transport_in_large_cities_worldwide.R
#
# The pipeline proceeds through 6 stages:
#   PART 0 - Setup (packages, working directory)
#   PART 1 - Data Inspection (JSON schema peek, field extraction)
#   PART 2 - JSON -> Parquet Conversion (safe reader handling NULLs, type coercion)
#   PART 3 - Network Reconstruction (per-city node & edge CSVs from parquet data)
#   PART 4 - Network Analytics (per-city metrics, summary table, .rds/.csv output)
#   PART 5 - Plotting (distribution plots, representative cities, cute plots)
#   PART 6 - Global Plots for LaTeX Appendix (one PDF page per metric)
#
# Outputs:
#   Data_parquet/              - Parquet files converted from raw JSON
#   city_graph_files/nodes/    - Per-city node CSVs
#   city_graph_files/edges/    - Per-city edge CSVs
#   city_graph_files/logs/     - Build summary log
#   analytics_out/             - per_city_summary.csv, citylines_metrics.rds
#   plots_out/                 - Plot PDFs (metric distributions, per-city)
#   cute_plots_out/            - Publication-quality n_nodes plots
#   global_plots/              - One PDF per metric for LaTeX appendix
#   json_peek_log.txt          - JSON schema inspection log
#   schema_dump.txt            - Extracted field schemas for all JSON files
# ==============================================================================


# ==============================================================================
# PART 0 - Setup: packages, working directory, reproducibility
# ==============================================================================

options(stringsAsFactors = FALSE)

# Set working directory to the script's own location when run via Rscript,
# so that relative paths (./Data, ./city_graph_files, etc.) work correctly.
.get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("^--file=", cmdArgs, value = TRUE)
  if (length(fileArg) == 0) return(getwd())
  normalizePath(sub("^--file=", "", fileArg), winslash = "/", mustWork = FALSE) |> dirname()
}
script_dir <- .get_script_dir()
setwd(script_dir)
message("Working directory set to: ", script_dir)

set.seed(1)

# Load all required packages (loaded once; no redundant library() calls below)
suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(arrow)
  library(sf)
  library(stringr)
  library(igraph)
  library(data.table)
  library(ggplot2)
  library(scales)
})


# ==============================================================================
# PART 1 - Data Inspection
#
# Preliminary overview of file sizes, JSON schema inspection, and field
# extraction. Writes json_peek_log.txt and schema_dump.txt for reference.
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 1: Data Inspection")
message(strrep("=", 70))

data_dir <- "./Data"
files <- list.files(data_dir, pattern = "\\.json$", full.names = TRUE)

# --- 1a. Overview of file sizes ---
file_overview <- tibble(
  file  = basename(files),
  path  = files,
  bytes = file.info(files)$size,
  mb    = round(bytes / 1024^2, 2)
) %>% arrange(desc(bytes))

cat("\nFile overview:\n")
print(file_overview)

# --- 1b. JSON schema peek (structure, top-level types, example values) ---
# Writes a human-readable log of each JSON file's structure.

peek_json <- function(path, n_show = 2) {
  x <- jsonlite::read_json(path, simplifyVector = FALSE)

  top <- if (is.list(x) && !is.null(names(x))) "object"
         else if (is.list(x)) "array"
         else class(x)[1]

  cat("\n============================\n")
  cat("File:", basename(path), "\n")
  cat("Path:", path, "\n")
  cat("Top-level:", top, "\n")

  if (top == "object") {
    cat("Keys:", paste(head(names(x), 50), collapse = ", "),
        if (length(names(x)) > 50) "..." else "", "\n")

    t <- map_chr(x, ~{
      if (is.null(.x)) "NULL"
      else if (is.atomic(.x)) paste0("atomic(", typeof(.x), ")")
      else if (is.list(.x) && !is.null(names(.x))) "object"
      else "array"
    })

    cat("\nFirst-level value types (counts):\n")
    print(head(sort(table(t), decreasing = TRUE), 20))

    cat("\nExample values (first keys):\n")
    for (i in seq_len(min(length(x), n_show))) {
      cat("\n--- Key:", names(x)[i], "---\n")
      str(x[[i]], max.level = 2)
    }

  } else if (top == "array") {
    cat("Length:", length(x), "\n")
    if (length(x) > 0) {
      cat("\nElement[1] structure:\n")
      str(x[[1]], max.level = 2)
      if (length(x) >= 2) {
        cat("\nElement[2] structure:\n")
        str(x[[2]], max.level = 2)
      }
    }
  } else {
    cat("Value:\n")
    print(x)
  }

  invisible(NULL)
}

write_peeks_to_log <- function(files, log_path = "json_peek_log.txt", n_show = 2) {
  cat("", file = log_path)
  cat("JSON peek log\n", file = log_path, append = TRUE)
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      file = log_path, append = TRUE)
  cat("Files:", length(files), "\n\n", file = log_path, append = TRUE)

  for (f in files) {
    out <- tryCatch(
      capture.output(peek_json(f, n_show = n_show)),
      error = function(e) {
        c("============================",
          paste("File:", basename(f)),
          paste("Path:", f),
          paste("ERROR:", conditionMessage(e)), "")
      }
    )
    cat(out, sep = "\n", file = log_path, append = TRUE)
    cat("\n", file = log_path, append = TRUE)
  }
  invisible(log_path)
}

log_file <- write_peeks_to_log(files, log_path = "json_peek_log.txt", n_show = 2)
cat("Wrote:", log_file, "\n")

# --- 1c. Quick schema + missing values + ranges (profiling) ---

profile_records <- function(path, sample_n = 2000) {
  x <- jsonlite::read_json(path, simplifyVector = FALSE)
  stopifnot(is.list(x), length(x) > 0)

  idx <- if (length(x) > sample_n) sample.int(length(x), sample_n) else seq_along(x)
  s <- x[idx]

  keys <- sort(unique(unlist(map(s, names))))
  cat("\nFile:", basename(path), "sampled:", length(s), "records\n")
  cat("Unique keys (sample):", length(keys), "\n")

  coverage <- tibble(
    key = keys,
    present_pct = map_dbl(keys, ~mean(map_lgl(s, \(r) .x %in% names(r)))) * 100
  ) %>% arrange(desc(present_pct))
  print(head(coverage, 30))

  type_summary <- map(keys, function(k) {
    vals <- map(s, \(r) r[[k]])
    tibble(
      key = k,
      type = names(sort(table(map_chr(vals, \(v)
        if (is.null(v)) "NULL"
        else if (is.atomic(v)) paste0("atomic(", typeof(v), ")")
        else if (is.list(v) && !is.null(names(v))) "object"
        else "array"
      )), decreasing = TRUE))[1]
    )
  }) %>% bind_rows()
  print(head(type_summary, 30))

  list(coverage = coverage, type_summary = type_summary)
}

profile_records(file.path(data_dir, "Stations.json"))

# --- 1d. Extract headers/fields for all JSON files -> schema_dump.txt ---

extract_schema <- function(path) {
  x <- jsonlite::read_json(path, simplifyVector = FALSE)

  title      <- if (!is.null(x$title)) as.character(x$title) else basename(path)
  fields     <- if (!is.null(x$fields)) unlist(x$fields) else character(0)
  type_names <- if (!is.null(x$type_names)) unlist(x$type_names) else NULL
  types      <- if (!is.null(x$types)) unlist(x$types) else NULL

  cat("\n============================\n")
  cat("File:", basename(path), "\n")
  cat("Title:", title, "\n")

  if (length(fields) == 0) {
    cat("No `fields` found. Top-level keys:\n")
    if (is.list(x) && !is.null(names(x))) cat(paste(names(x), collapse = ", "), "\n")
    return(invisible(NULL))
  }

  cat("n_fields:", length(fields), "\n")
  cat("fields:\n")
  cat(paste0(" - ", fields), sep = "\n")
  cat("\n")

  if (!is.null(type_names) && length(type_names) == length(fields)) {
    cat("type_names:\n")
    for (i in seq_along(fields)) {
      cat(sprintf(" - %-25s %s\n", fields[i], type_names[i]))
    }
  } else if (!is.null(types) && length(types) == length(fields)) {
    cat("types:\n")
    for (i in seq_along(fields)) {
      cat(sprintf(" - %-25s %s\n", fields[i], types[i]))
    }
  }

  invisible(NULL)
}

json_files <- list.files(data_dir, pattern = "\\.json$", full.names = TRUE)
sink("schema_dump.txt")
walk(json_files, extract_schema)
sink()
cat("Wrote schema_dump.txt\n")


# ==============================================================================
# PART 2 - JSON -> Parquet Conversion
#
# The JSON files are basically tabular (no deep/irregular nesting), so we
# convert them to Parquet for faster loading and better memory efficiency.
# Some JSON rows contain NULL cells which break bind_rows(); we replace
# NULLs with NA to handle this robustly.
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 2: JSON -> Parquet Conversion")
message(strrep("=", 70))

# Replace NULLs inside a list with NA (keeps list length stable for bind_rows)
null_to_na <- function(v) {
  map(v, ~ if (is.null(.x)) NA else .x)
}

# Safe JSON reader: handles NULL values that would otherwise break bind_rows()
read_columnar_json_safe <- function(path) {
  x <- jsonlite::read_json(path, simplifyVector = FALSE)
  stopifnot(is.list(x), !is.null(x$fields), !is.null(x$values))

  fields <- unlist(x$fields)
  p <- length(fields)

  rows <- map(x$values, function(r) {
    if (is.null(r)) r <- vector("list", p)
    if (!is.list(r)) r <- as.list(r)

    # Convert NULL cells to NA
    r <- null_to_na(r)

    # Pad/truncate to match fields length
    if (length(r) < p) r[(length(r) + 1):p] <- NA
    if (length(r) > p) r <- r[seq_len(p)]

    names(r) <- fields
    r
  })

  df <- bind_rows(rows) %>% as_tibble()

  attr(df, "json_meta") <- list(
    title      = x$title,
    types      = x$types,
    type_names = x$type_names
  )

  df
}

# Type coercion based on type_names metadata embedded in the JSON
coerce_from_type_names <- function(df) {
  meta <- attr(df, "json_meta")
  if (is.null(meta) || is.null(meta$type_names)) return(df)

  type_names <- unlist(meta$type_names)
  fields <- names(df)

  n <- min(length(fields), length(type_names))
  fields <- fields[seq_len(n)]
  type_names <- tolower(type_names[seq_len(n)])

  for (i in seq_len(n)) {
    col <- fields[i]
    tn  <- type_names[i]

    if (tn %in% c("integer", "int", "int4", "int8", "bigint", "smallint")) {
      df[[col]] <- suppressWarnings(as.integer(df[[col]]))
    } else if (tn %in% c("double", "float", "float4", "float8", "numeric", "real", "decimal")) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    } else if (tn %in% c("boolean", "bool")) {
      df[[col]] <- as.logical(df[[col]])
    } else if (tn %in% c("text", "string", "character", "varchar")) {
      df[[col]] <- as.character(df[[col]])
    }
    # Leave unknown types (e.g., WKT geometry) as-is
  }

  df
}

# Convert all JSON files in Data/ to Parquet in Data_parquet/
out_dir <- "./Data_parquet"
dir.create(out_dir, showWarnings = FALSE)

json_files <- list.files(data_dir, pattern = "\\.json$", full.names = TRUE)

purrr::walk(json_files, function(f) {
  message("Converting: ", basename(f))
  df <- read_columnar_json_safe(f)
  df <- coerce_from_type_names(df)
  out_path <- file.path(out_dir, sub("\\.json$", ".parquet", basename(f)))
  arrow::write_parquet(df, out_path)
})

message("Done. Wrote ", length(json_files), " parquet files to: ", out_dir)


# ==============================================================================
# PART 3 - Network Reconstruction (per-city node & edge CSVs)
#
# Build per-city node and edge lists from the public-transport dataset.
#
# Nodes: nodeID, nodeLabel, latitude, longitude, mode, year
# Edges: nodeID_from, nodeID_to, mode, line, year
#
# Conventions:
# - Node IDs remapped 1..N_city per city
# - Edges are (section, line) => multigraph
# - Edge line = Lines.url_name
# - Mode for nodes = concatenation of supported modes with "|"
# - Year = opening year
#     Nodes: Stations.opening, fallback to earliest station_lines.fromyear
#     Edges: section_lines.fromyear, fallback to Sections.opening
# - Section endpoints snapped to nearest stations using great-circle distance
# - Per-city snap threshold = max(200m, q99 of distances), capped at 1500m
# - Self-loops dropped
# - Undirected edges: output nodeID_from = min, nodeID_to = max
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 3: Network Reconstruction")
message(strrep("=", 70))

# --- Config ---
PARQ_DIR  <- "./Data_parquet"
GRAPH_DIR <- "./city_graph_files"

THR_MIN_M <- 200     # minimum snap threshold (metres)
THR_Q     <- 0.99    # quantile of endpoint distances for adaptive threshold
THR_CAP_M <- 1500    # maximum snap threshold cap (metres)

CRS_WGS84 <- 4326    # EPSG code for WGS84 (lon/lat)

dir.create(GRAPH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(GRAPH_DIR, "nodes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(GRAPH_DIR, "edges"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(GRAPH_DIR, "logs"),  showWarnings = FALSE, recursive = TRUE)

# --- Helpers for reconstruction ---

safe_city_slug <- function(x) {
  x <- ifelse(is.na(x) | x == "", "unknown", x)
  x <- str_replace_all(x, "[^A-Za-z0-9_-]+", "-")
  x <- str_replace_all(x, "-+", "-")
  x <- str_replace_all(x, "^-|-$", "")
  tolower(x)
}

to_int <- function(x) suppressWarnings(as.integer(x))

# Extract lon/lat from "POINT(lon lat)" WKT strings
parse_point_wkt_lonlat <- function(wkt_vec) {
  wkt_vec <- as.character(wkt_vec)
  m <- str_match(wkt_vec, "POINT\\s*\\(\\s*([-0-9\\.eE]+)\\s+([-0-9\\.eE]+)\\s*\\)")
  data.frame(
    lon = suppressWarnings(as.numeric(m[, 2])),
    lat = suppressWarnings(as.numeric(m[, 3]))
  )
}

make_sf_points <- function(lon, lat) {
  st_as_sf(
    data.frame(lon = lon, lat = lat),
    coords = c("lon", "lat"),
    crs = CRS_WGS84,
    remove = FALSE
  )
}

# Extract first and last coordinates of each section geometry
section_endpoints_sf <- function(sections_df) {
  geom <- st_as_sfc(sections_df$geometry, crs = CRS_WGS84)
  coords_list <- lapply(seq_along(geom), function(i) {
    xy <- st_coordinates(geom[i])
    if (nrow(xy) == 0) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    c(xy[1, "X"], xy[1, "Y"], xy[nrow(xy), "X"], xy[nrow(xy), "Y"])
  })
  mat <- do.call(rbind, coords_list)
  out <- sections_df %>%
    transmute(
      section_id = id,
      x_from = mat[, 1], y_from = mat[, 2],
      x_to   = mat[, 3], y_to   = mat[, 4]
    )
  out
}

# Snap section endpoints to nearest stations within a threshold
snap_endpoints_to_stations <- function(endpoints_df, stations_sf, thr_m) {
  pts_from <- make_sf_points(endpoints_df$x_from, endpoints_df$y_from)
  pts_to   <- make_sf_points(endpoints_df$x_to,   endpoints_df$y_to)

  idx_from <- st_nearest_feature(pts_from, stations_sf)
  idx_to   <- st_nearest_feature(pts_to,   stations_sf)

  d_from <- as.numeric(st_distance(pts_from, stations_sf[idx_from, ], by_element = TRUE))
  d_to   <- as.numeric(st_distance(pts_to,   stations_sf[idx_to, ],   by_element = TRUE))

  from_station_id <- stations_sf$station_id[idx_from]
  to_station_id   <- stations_sf$station_id[idx_to]

  endpoints_df %>%
    mutate(
      from_station_id = from_station_id,
      to_station_id   = to_station_id,
      d_from_m = d_from,
      d_to_m   = d_to,
      ok = !is.na(d_from_m) & !is.na(d_to_m) &
           d_from_m <= thr_m & d_to_m <= thr_m &
           !is.na(from_station_id) & !is.na(to_station_id) &
           (from_station_id != to_station_id)
    )
}

# --- Load parquet files ---
read_p <- function(name) read_parquet(file.path(PARQ_DIR, paste0(name, ".parquet")))

cities          <- read_p("Cities")
lines           <- read_p("Lines")
sections        <- read_p("Sections")
stations        <- read_p("Stations")
systems         <- read_p("Systems")
transport_modes <- read_p("Transport_modes")
section_lines   <- read_p("Section_lines")
station_lines   <- read_p("Station_lines")

# Ensure key columns are integer
cities$id                  <- to_int(cities$id)
lines$id                   <- to_int(lines$id)
lines$city_id              <- to_int(lines$city_id)
lines$system_id            <- to_int(lines$system_id)
lines$transport_mode_id    <- to_int(lines$transport_mode_id)

stations$id                <- to_int(stations$id)
stations$city_id           <- to_int(stations$city_id)
stations$opening           <- to_int(stations$opening)

sections$id                <- to_int(sections$id)
sections$city_id           <- to_int(sections$city_id)
sections$opening           <- to_int(sections$opening)

section_lines$section_id   <- to_int(section_lines$section_id)
section_lines$line_id      <- to_int(section_lines$line_id)
section_lines$city_id      <- to_int(section_lines$city_id)
section_lines$fromyear     <- to_int(section_lines$fromyear)

station_lines$station_id   <- to_int(station_lines$station_id)
station_lines$line_id      <- to_int(station_lines$line_id)
station_lines$city_id      <- to_int(station_lines$city_id)
station_lines$fromyear     <- to_int(station_lines$fromyear)

transport_modes$id         <- to_int(transport_modes$id)

# Pre-join lines -> modes (for speed in the per-city loop)
lines_modes <- lines %>%
  left_join(transport_modes %>% select(transport_mode_id = id, mode_name = name),
            by = "transport_mode_id") %>%
  mutate(
    line_key  = ifelse(!is.na(url_name) & url_name != "", url_name, name),
    line_key  = as.character(line_key),
    mode_name = as.character(mode_name)
  ) %>%
  select(line_id = id, city_id, line_key, mode_name)

# --- Main per-city loop: build node & edge CSVs ---
log_rows <- list()
all_city_ids <- sort(unique(cities$id))

for (cid in all_city_ids) {
  tryCatch({

    city_row   <- cities %>% filter(id == cid) %>% slice(1)
    city_slug  <- safe_city_slug(city_row$url_name %||% city_row$name %||% as.character(cid))
    city_label <- city_row$name %||% city_slug

    message("\n=== City ", cid, " (", city_label, ") ===")

    # Filter city-specific tables
    st_city  <- stations %>% filter(city_id == cid)
    if (nrow(st_city) == 0) {
      message("No stations; skipping city.")
      next
    }

    se_city  <- sections %>% filter(city_id == cid)
    sl_city  <- section_lines %>% filter(city_id == cid)
    stl_city <- station_lines %>% filter(city_id == cid)

    # Build station lon/lat
    lonlat <- parse_point_wkt_lonlat(st_city$geometry)
    st_city <- st_city %>%
      mutate(lon = lonlat$lon, lat = lonlat$lat)

    # Node year: Stations.opening, else earliest station_lines.fromyear
    st_year_fallback <- stl_city %>%
      filter(!is.na(fromyear)) %>%
      group_by(station_id) %>%
      summarise(year_fallback = suppressWarnings(min(fromyear, na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(year_fallback = ifelse(is.infinite(year_fallback),
                                    NA_integer_, as.integer(year_fallback)))

    # Node mode: concatenate unique modes served by station (via station_lines -> lines -> modes)
    st_modes <- stl_city %>%
      inner_join(lines_modes %>% filter(city_id == cid), by = c("line_id", "city_id")) %>%
      filter(!is.na(mode_name) & mode_name != "") %>%
      group_by(station_id) %>%
      summarise(mode = paste(sort(unique(mode_name)), collapse = "|"), .groups = "drop")

    st_city <- st_city %>%
      left_join(st_year_fallback, by = c("id" = "station_id")) %>%
      left_join(st_modes, by = c("id" = "station_id")) %>%
      mutate(
        year      = ifelse(!is.na(opening), opening, year_fallback),
        nodeLabel = ifelse(!is.na(name) & name != "", as.character(name), as.character(id)),
        mode      = ifelse(is.na(mode) | mode == "", NA_character_, mode)
      )

    # Remap node IDs 1..N_city
    st_city <- st_city %>%
      arrange(id) %>%
      mutate(nodeID = row_number())

    station_id_to_nodeID <- st_city %>%
      select(station_id = id, nodeID)

    # Write nodes CSV
    nodes_out <- st_city %>%
      select(nodeID, nodeLabel, latitude = lat, longitude = lon, mode, year)

    nodes_path <- file.path(GRAPH_DIR, "nodes",
                            paste0("nodes_city_", cid, "_", city_slug, ".csv"))
    write.csv(nodes_out, nodes_path, row.names = FALSE, na = "")

    # --- Build section endpoints and snap to stations ---
    if (nrow(se_city) == 0 || nrow(sl_city) == 0) {
      message("No sections or no section_lines; writing nodes only.")
      log_rows[[length(log_rows) + 1]] <- data.frame(
        city_id = cid, city = city_label,
        n_nodes = nrow(nodes_out), n_sections = nrow(se_city),
        n_section_lines = nrow(sl_city), threshold_m = NA_real_,
        edges_written = 0L, edges_dropped_unmatched = nrow(sl_city),
        edges_dropped_selfloops = 0L, stringsAsFactors = FALSE
      )
      next
    }

    # Build sf stations for snapping
    stations_sf <- st_city %>%
      mutate(lon = as.numeric(lon), lat = as.numeric(lat)) %>%
      filter(is.finite(lon), is.finite(lat)) %>%
      transmute(station_id = id, lon, lat) %>%
      st_as_sf(coords = c("lon", "lat"), crs = CRS_WGS84, remove = FALSE)

    if (nrow(stations_sf) == 0) {
      message("No valid station coordinates after parsing; writing nodes only.")
      log_rows[[length(log_rows) + 1]] <- data.frame(
        city_id = cid, city = city_label,
        n_nodes = nrow(nodes_out), n_sections = nrow(se_city),
        n_section_lines = nrow(sl_city), threshold_m = NA_real_,
        edges_written = 0L, edges_dropped_unmatched = nrow(sl_city),
        edges_dropped_selfloops = 0L, stringsAsFactors = FALSE
      )
      next
    }

    # Extract section endpoints
    endpoints <- section_endpoints_sf(se_city)

    # Compute preliminary nearest distances to derive adaptive threshold
    tmp_snapped <- snap_endpoints_to_stations(endpoints, stations_sf, thr_m = Inf)

    d_all <- c(tmp_snapped$d_from_m, tmp_snapped$d_to_m)
    d_all <- d_all[is.finite(d_all) & !is.na(d_all)]

    if (length(d_all) == 0) {
      message("Could not compute endpoint distances; skipping edges for this city.")
      log_rows[[length(log_rows) + 1]] <- data.frame(
        city_id = cid, city = city_label,
        n_nodes = nrow(nodes_out), n_sections = nrow(se_city),
        n_section_lines = nrow(sl_city), threshold_m = NA_real_,
        edges_written = 0L, edges_dropped_unmatched = nrow(sl_city),
        edges_dropped_selfloops = 0L, stringsAsFactors = FALSE
      )
      next
    }

    thr_city <- as.numeric(quantile(d_all, probs = THR_Q, na.rm = TRUE))
    thr_city <- max(THR_MIN_M, thr_city)
    thr_city <- min(THR_CAP_M, thr_city)

    message("Threshold (m): ", round(thr_city, 2),
            " (q", THR_Q * 100, ", min ", THR_MIN_M, ", cap ", THR_CAP_M, ")")

    snapped <- snap_endpoints_to_stations(endpoints, stations_sf, thr_m = thr_city)

    # Keep only valid snaps (also drops self-loops)
    snapped_ok <- snapped %>%
      filter(ok) %>%
      select(section_id, from_station_id, to_station_id)

    # Build edges as (section, line) pairs
    edges0 <- sl_city %>%
      select(section_id, line_id, fromyear) %>%
      inner_join(lines_modes %>% filter(city_id == cid), by = c("line_id")) %>%
      left_join(se_city %>% select(section_id = id, section_opening = opening),
                by = "section_id") %>%
      mutate(
        year = ifelse(!is.na(fromyear), fromyear, section_opening),
        line = line_key,
        mode = mode_name
      ) %>%
      select(section_id, line, mode, year)

    # Attach snapped endpoints
    edges1 <- edges0 %>%
      left_join(snapped_ok, by = "section_id") %>%
      inner_join(station_id_to_nodeID, by = c("from_station_id" = "station_id")) %>%
      rename(nodeID_from = nodeID) %>%
      inner_join(station_id_to_nodeID, by = c("to_station_id" = "station_id")) %>%
      rename(nodeID_to = nodeID)

    # Undirected canonical ordering + drop any remaining self-loops
    edges2 <- edges1 %>%
      mutate(
        nodeID_min = pmin(nodeID_from, nodeID_to),
        nodeID_max = pmax(nodeID_from, nodeID_to)
      ) %>%
      filter(nodeID_min != nodeID_max) %>%
      transmute(
        nodeID_from = nodeID_min,
        nodeID_to   = nodeID_max,
        mode, line, year
      )

    edges_path <- file.path(GRAPH_DIR, "edges",
                            paste0("edges_city_", cid, "_", city_slug, ".csv"))
    write.csv(edges2, edges_path, row.names = FALSE, na = "")

    # Logging
    matched_sections <- unique(snapped_ok$section_id)
    n_unmatched_section_lines <- sum(!edges0$section_id %in% matched_sections)
    n_selfloop_sections <- sum(tmp_snapped$from_station_id == tmp_snapped$to_station_id,
                               na.rm = TRUE)

    log_rows[[length(log_rows) + 1]] <- data.frame(
      city_id = cid, city = city_label,
      n_nodes = nrow(nodes_out), n_sections = nrow(se_city),
      n_section_lines = nrow(sl_city), threshold_m = thr_city,
      edges_written = nrow(edges2),
      edges_dropped_unmatched = n_unmatched_section_lines,
      edges_dropped_selfloops = n_selfloop_sections,
      stringsAsFactors = FALSE
    )

    message("Wrote: ", basename(nodes_path), " (", nrow(nodes_out), " nodes)")
    message("Wrote: ", basename(edges_path), " (", nrow(edges2), " edges)")

  }, error = function(e) {
    message("City ", cid, " FAILED: ", conditionMessage(e))
    log_rows[[length(log_rows) + 1]] <<- data.frame(
      city_id = cid, city = NA_character_,
      n_nodes = NA_integer_, n_sections = NA_integer_,
      n_section_lines = NA_integer_, threshold_m = NA_real_,
      edges_written = NA_integer_, edges_dropped_unmatched = NA_integer_,
      edges_dropped_selfloops = NA_integer_, stringsAsFactors = FALSE
    )
  })
}

# Write reconstruction summary log
log_df   <- bind_rows(log_rows)
log_path <- file.path(GRAPH_DIR, "logs", "build_summary.csv")
write.csv(log_df, log_path, row.names = FALSE)
message("Output directory: ", GRAPH_DIR)
message("Summary log: ", log_path)


# ==============================================================================
# PART 4 - Network Analytics
#
# Build per-city network summaries from nodes_*.csv + edges_*.csv
# (multigraph + simplified graph), then compute over-cities distributions.
# Saves .rds and .csv to analytics_out/.
#
# Reuses key patterns from combined.R (Prof. Manlio De Domenico):
# - normalize_vec, vec2pal
# - read.csv -> edges df -> graph_from_data_frame -> simplify
# - lon/lat layout matrix
# - self-loop filtering
# - community detection + modularity
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 4: Network Analytics")
message(strrep("=", 70))

# --- Config ---
NODES_DIR <- "city_graph_files/nodes"
EDGES_DIR <- "city_graph_files/edges"
ANALYTICS_OUT <- "analytics_out"

PATH_SAMPLE_N        <- 400   # nodes sampled for distance computations
PATH_MAX_COMPONENT_N <- 5000  # if giant component larger, sample within it

dir.create(ANALYTICS_OUT, showWarnings = FALSE, recursive = TRUE)

# --- Reused helpers from common.R (Prof. Manlio De Domenico) ---

normalize_vec <- function(v_centr) {
  (v_centr - min(v_centr, na.rm = TRUE)) /
    (max(v_centr, na.rm = TRUE) - min(v_centr, na.rm = TRUE))
}

vec2pal <- function(v_centr, mypal) {
  v <- v_centr
  v[is.na(v)] <- min(v, na.rm = TRUE)
  val_idxs <- 1 + floor((length(mypal) - 1) * normalize_vec(v))
  val_idxs[val_idxs < 1] <- 1
  val_idxs[val_idxs > length(mypal)] <- length(mypal)
  mypal[val_idxs]
}

# --- Citylines-specific helpers ---

# Parse "mode" field which is "a|b|c" (pipe-separated, alphabetical)
split_modes <- function(mode_str) {
  if (is.na(mode_str) || nchar(mode_str) == 0) return(character(0))
  m <- unlist(strsplit(tolower(as.character(mode_str)), "\\|"))
  m <- trimws(m)
  m[m != ""]
}

# Build multigraph and simplified graph from node + edge data.tables
build_graphs <- function(nodes_dt, edges_dt) {
  # Remove self-loops early
  edges_dt <- edges_dt[nodeID_from != nodeID_to]

  vdf <- data.frame(name = as.character(nodes_dt$nodeID), stringsAsFactors = FALSE)

  g_multi <- igraph::graph_from_data_frame(
    d = data.frame(from = as.character(edges_dt$nodeID_from),
                   to   = as.character(edges_dt$nodeID_to),
                   stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = vdf
  )

  # Simple version for shortest paths / communities
  g_simple <- igraph::simplify(g_multi, remove.multiple = TRUE, remove.loops = TRUE)

  list(g_multi = g_multi, g_simple = g_simple, edges_clean = edges_dt)
}

# Layout matrix from lon/lat (for geo-spatial network plots)
layout_from_nodes <- function(nodes_dt, g) {
  idx <- match(V(g)$name, as.character(nodes_dt$nodeID))
  cbind(nodes_dt$longitude[idx], nodes_dt$latitude[idx])
}

# Compute shortest-path stats on the giant component (sampled for large graphs)
path_stats <- function(g_simple) {
  if (vcount(g_simple) < 2 || ecount(g_simple) < 1) {
    return(list(mean_dist = NA_real_, median_dist = NA_real_,
                p90_dist = NA_real_, diameter = NA_real_))
  }

  comps    <- components(g_simple)
  giant_id <- which.max(comps$csize)
  v_giant  <- V(g_simple)[comps$membership == giant_id]
  g_giant  <- induced_subgraph(g_simple, v_giant)

  nG <- vcount(g_giant)
  if (nG < 2) {
    return(list(mean_dist = NA_real_, median_dist = NA_real_,
                p90_dist = NA_real_, diameter = NA_real_))
  }

  # Sample nodes for distance computations
  sample_n <- min(PATH_SAMPLE_N, nG)
  if (nG > PATH_MAX_COMPONENT_N) {
    vs <- sample(V(g_giant), sample_n)
  } else {
    vs <- if (sample_n < nG) sample(V(g_giant), sample_n) else V(g_giant)
  }

  D    <- distances(g_giant, v = vs, to = V(g_giant))
  dvec <- as.numeric(D)
  dvec <- dvec[is.finite(dvec) & dvec > 0]

  if (length(dvec) == 0) {
    return(list(mean_dist = NA_real_, median_dist = NA_real_,
                p90_dist = NA_real_, diameter = NA_real_))
  }

  diam <- tryCatch(diameter(g_giant, directed = FALSE, unconnected = TRUE),
                   error = function(e) NA_real_)

  list(
    mean_dist   = mean(dvec),
    median_dist = median(dvec),
    p90_dist    = as.numeric(quantile(dvec, 0.90)),
    diameter    = diam
  )
}

# Community detection (mesoscale structure)
community_stats <- function(g_simple) {
  if (vcount(g_simple) < 2 || ecount(g_simple) < 1) {
    return(list(mod_louvain = NA_real_, k_louvain = NA_integer_,
                mod_infomap = NA_real_, k_infomap = NA_integer_))
  }

  lou <- tryCatch(cluster_louvain(g_simple, resolution = 1), error = function(e) NULL)
  inf <- tryCatch(cluster_infomap(g_simple), error = function(e) NULL)

  mod_lou <- if (!is.null(lou)) modularity(lou)                   else NA_real_
  k_lou   <- if (!is.null(lou)) length(unique(membership(lou)))   else NA_integer_
  mod_inf <- if (!is.null(inf)) modularity(inf)                   else NA_real_
  k_inf   <- if (!is.null(inf)) length(unique(membership(inf)))   else NA_integer_

  list(mod_louvain = mod_lou, k_louvain = k_lou,
       mod_infomap = mod_inf, k_infomap = k_inf)
}

# --- Main: loop over cities, compute summary table ---

node_files <- list.files(NODES_DIR, pattern = "^nodes_.*\\.csv$", full.names = TRUE)
if (length(node_files) == 0) stop("No node files found in ", NODES_DIR)

summaries <- list()

for (nf in node_files) {
  base     <- basename(nf)
  city_tag <- sub("^nodes_city_", "", sub("\\.csv$", "", base))

  # Parse city_index and city_name from filename
  m <- str_match(city_tag, "^([0-9]+)_(.*)$")
  city_index <- if (!is.na(m[1, 2])) as.integer(m[1, 2]) else NA_integer_
  city_name  <- if (!is.na(m[1, 3])) m[1, 3]             else city_tag

  ef <- file.path(EDGES_DIR, paste0("edges_city_", city_tag, ".csv"))
  if (!file.exists(ef)) {
    message("Skipping (missing edges): ", city_tag)
    next
  }

  nodes <- fread(nf)
  edges <- fread(ef)

  # Basic validation
  required_nodes <- c("nodeID", "nodeLabel", "latitude", "longitude", "mode", "year")
  required_edges <- c("nodeID_from", "nodeID_to", "mode", "line", "year")
  if (!all(required_nodes %in% names(nodes))) stop("Bad node headers in ", nf)
  if (!all(required_edges %in% names(edges))) stop("Bad edge headers in ", ef)

  # Build graphs
  gs       <- build_graphs(nodes, edges)
  g_multi  <- gs$g_multi
  g_simple <- gs$g_simple

  # Basic metrics
  n        <- vcount(g_simple)
  m_multi  <- ecount(g_multi)
  m_simple <- ecount(g_simple)

  dens           <- if (n > 1) edge_density(g_simple, loops = FALSE) else NA_real_
  avg_deg_simple <- if (n > 0) mean(degree(g_simple))                else NA_real_
  avg_deg_multi  <- if (vcount(g_multi) > 0) mean(degree(g_multi))   else NA_real_

  comps      <- components(g_simple)
  n_comp     <- comps$no
  giant_frac <- if (n > 0) max(comps$csize) / n else NA_real_

  trans <- tryCatch(transitivity(g_simple, type = "globalundirected"),
                    error = function(e) NA_real_)

  # Path statistics (sampled)
  ps <- path_stats(g_simple)

  # Community structure (mesoscale)
  cs <- community_stats(g_simple)

  # Multigraph multiplicity: how many parallel edges per unique pair?
  el <- as.data.table(as_edgelist(g_multi, names = TRUE))
  if (nrow(el) == 0) {
    mean_mult <- NA_real_
    max_mult  <- NA_real_
  } else {
    setnames(el, c("V1", "V2"), c("u", "v"))
    el[, `:=`(a = pmin(u, v), b = pmax(u, v))]
    mult_counts <- el[, .N, by = .(a, b)]$N
    mean_mult   <- mean(mult_counts)
    max_mult    <- max(mult_counts)
  }

  summaries[[city_tag]] <- data.table(
    city_tag    = city_tag,
    city_index  = city_index,
    city_name   = city_name,

    n_nodes     = n,
    m_edges_multi  = m_multi,
    m_edges_simple = m_simple,

    density_simple     = dens,
    avg_degree_simple  = avg_deg_simple,
    avg_degree_multi   = avg_deg_multi,

    n_components          = n_comp,
    giant_component_frac  = giant_frac,
    transitivity          = trans,

    mean_distance   = ps$mean_dist,
    median_distance = ps$median_dist,
    p90_distance    = ps$p90_dist,
    diameter        = ps$diameter,

    modularity_louvain     = cs$mod_louvain,
    n_communities_louvain  = cs$k_louvain,
    modularity_infomap     = cs$mod_infomap,
    n_communities_infomap  = cs$k_infomap,

    mean_edge_multiplicity = mean_mult,
    max_edge_multiplicity  = max_mult
  )
}

summary_dt <- rbindlist(summaries, fill = TRUE)
setorder(summary_dt, city_index, city_name)

# Save summary outputs
fwrite(summary_dt, file.path(ANALYTICS_OUT, "per_city_summary.csv"))

results_obj <- list(
  per_city_summary = summary_dt,
  config = list(
    nodes_dir     = NODES_DIR,
    edges_dir     = EDGES_DIR,
    path_sample_n = PATH_SAMPLE_N
  ),
  plot_functions = c(
    "plot_degree_distribution(g_multi, g_simple, city_tag)",
    "plot_network_on_map(g_simple, nodes_dt, city_tag, color_by)",
    "plot_metric_distribution(summary_dt, metric)",
    "pick_representative_cities(summary_dt, metric)"
  )
)

saveRDS(results_obj, file.path(ANALYTICS_OUT, "citylines_metrics.rds"))

message("Done. Wrote: ",
        file.path(ANALYTICS_OUT, "per_city_summary.csv"),
        " and ",
        file.path(ANALYTICS_OUT, "citylines_metrics.rds"))


# ==============================================================================
# PART 5 - Plotting
#
# Plotting utilities (for both global and per-city individual plots):
# - Degree distribution (multigraph vs. simple)
# - Network-on-map (colored by Louvain communities, sized by degree)
# - Over-cities metric distribution histograms
# - Representative city selection (closest to q10/mean/q90)
# Then:
# - Quick single-metric plot
# - Full loop over all metrics (filtered cities with avg degree >= 1)
# - Cute plots: n_nodes distribution + reference city graphs (for paper)
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 5: Plotting")
message(strrep("=", 70))

# --- Plotting utility functions ---

plot_degree_distribution <- function(g_multi, g_simple, city_tag) {
  d_multi  <- degree(g_multi)
  d_simple <- degree(g_simple)

  dt <- rbind(
    data.table(degree = d_multi,  type = "multigraph"),
    data.table(degree = d_simple, type = "simple")
  )

  ggplot(dt, aes(x = degree)) +
    geom_histogram(bins = 40) +
    facet_wrap(~type, nrow = 1) +
    labs(title = paste0("Degree distribution - ", city_tag),
         x = "degree", y = "count") +
    theme_minimal()
}

plot_network_on_map <- function(g_simple, nodes_dt, city_tag, resolution = 1) {
  lay <- layout_from_nodes(nodes_dt, g_simple)

  # Size by degree
  deg   <- degree(g_simple)
  vsize <- 1 + 4 * (deg / max(deg))

  # Color by community (Louvain)
  comm <- cluster_louvain(g_simple, resolution = resolution)
  memb <- membership(comm)
  pal  <- hcl.colors(length(unique(memb)), "Set3")
  cols <- pal[as.integer(as.factor(memb))]

  plot(
    g_simple, layout = lay,
    vertex.size  = vsize,
    vertex.label = NA,
    vertex.color = cols,
    edge.color   = alpha("black", 0.15),
    main         = paste0(city_tag, " (community color, degree size)")
  )
}

# Over-cities distribution histogram
plot_metric_distribution <- function(summary_dt, metric) {
  stopifnot(metric %in% names(summary_dt))
  ggplot(summary_dt[!is.na(get(metric))], aes(x = get(metric))) +
    geom_histogram(bins = 35) +
    labs(title = paste0("Across-cities distribution: ", metric),
         x = metric, y = "cities") +
    theme_minimal()
}

# Robust numeric mode (handles ties)
.numeric_mode <- function(x, ties = c("first", "min", "max", "mean")) {
  ties <- match.arg(ties)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)

  tab   <- table(x)
  max_n <- max(tab)
  modes <- as.numeric(names(tab)[tab == max_n])

  if (length(modes) == 1) return(modes)

  switch(ties,
    first = modes[1],
    min   = min(modes),
    max   = max(modes),
    mean  = mean(modes)
  )
}

# Pick representative city closest to a reference value for a given metric.
# ref can be:
#   - numeric in (0,1): percentile (e.g., 0.1, 0.9)
#   - "mean" or "mode"
pick_representative_cities <- function(summary_dt, metric,
                                        ref = c("mean", 0.25, 0.75),
                                        mode_ties = "first") {
  dt <- summary_dt[!is.na(get(metric))]
  if (nrow(dt) < 1) return(dt[0])

  x <- dt[[metric]]

  target    <- NULL
  ref_label <- NULL

  if (is.numeric(ref) && length(ref) == 1) {
    if (!(ref > 0 && ref < 1))
      stop("If numeric, `ref` must be a single value in (0, 1), e.g. 0.1.")
    target    <- as.numeric(quantile(x, probs = ref, na.rm = TRUE, names = FALSE))
    ref_label <- paste0("p", round(ref * 100))
  } else if (is.character(ref) && length(ref) == 1) {
    r <- tolower(ref)
    if (r == "mean") {
      target    <- mean(x, na.rm = TRUE)
      ref_label <- "mean"
    } else if (r == "mode") {
      target    <- .numeric_mode(x, ties = mode_ties)
      ref_label <- "mode"
    } else {
      stop("Unknown ref='", ref, "'. Use numeric percentile (0,1), 'mean', or 'mode'.")
    }
  } else {
    stop("`ref` must be a single numeric percentile in (0,1) or 'mean'/'mode'.")
  }

  if (is.na(target)) return(dt[0])

  dt[, dist_ref := abs(get(metric) - target)]
  dt[which.min(dist_ref),
     .(city_tag, city_index, city_name,
       metric_value = get(metric),
       ref = ref_label,
       ref_value = target)]
}

# --- 5a. Quick single-metric example plot ---

dir.create("plots_out", showWarnings = FALSE, recursive = TRUE)

pdf("plots_out/vscode_plots.pdf", width = 9, height = 5, onefile = TRUE)

res <- readRDS("analytics_out/citylines_metrics.rds")
S   <- res$per_city_summary

# Distribution over cities for density_simple
p <- plot_metric_distribution(S, "density_simple"); print(p)

# Representative cities (mean, Q1, Q3) for that metric
reps <- pick_representative_cities(S, "density_simple"); print(reps)

# Per-city degree distribution for the first representative city
city_tag <- reps$city_tag[1]
nodes <- fread(file.path(NODES_DIR, paste0("nodes_city_", city_tag, ".csv")))
edges <- fread(file.path(EDGES_DIR, paste0("edges_city_", city_tag, ".csv")))
gs <- build_graphs(nodes, edges)
print(plot_degree_distribution(gs$g_multi, gs$g_simple, city_tag))

dev.off()

# --- 5b. Full loop over all metrics (filtered: avg degree >= 1) ---

pdf("plots_out/vscode_plots_filtered_avgdeggt1.pdf", width = 9, height = 5, onefile = TRUE)

res <- readRDS("analytics_out/citylines_metrics.rds")
S   <- res$per_city_summary

unfilter_n_cities <- nrow(S)

# Filter cities with avg degree >= 1
deg_col <- if ("avg_degree_simple" %in% names(S)) {
  "avg_degree_simple"
} else if ("avg_degree_multi" %in% names(S)) {
  "avg_degree_multi"
} else {
  stop("No avg_degree_simple or avg_degree_multi column found in per_city_summary.")
}

S <- S[!is.na(get(deg_col)) & get(deg_col) >= 1.0]
filtered_n_cities <- nrow(S)
message("Cities kept after degree filter (", deg_col, " >= 1): ",
        filtered_n_cities, " (", round(filtered_n_cities / unfilter_n_cities * 100), "%)")

if (nrow(S) == 0) stop("No cities left after filtering by average degree >= 1.")

# Identify numeric metric columns (exclude identifiers)
id_cols     <- c("city_tag", "city_index", "city_name")
metric_cols <- setdiff(names(S), id_cols)
metric_cols <- metric_cols[sapply(S[, ..metric_cols], is.numeric)]

# Drop metrics that are all NA or degenerate (< 5 values or < 2 unique)
metric_cols <- metric_cols[
  sapply(metric_cols, function(m) {
    x <- S[[m]]
    x <- x[!is.na(x)]
    length(x) >= 5 && length(unique(x)) >= 2
  })
]

message("Will plot metrics: ", paste(metric_cols, collapse = ", "))

# Loop over all metrics
for (metric in metric_cols) {
  cat("\n========================\n")
  cat("Metric:", metric, "\n")
  cat("========================\n")

  # 1) Distribution over cities (filtered set)
  p <- plot_metric_distribution(S, metric)
  print(p)

  # 2) Representative cities (q10, mean, q90)
  reps <- rbind(
    pick_representative_cities(S, metric, ref = 0.1)[, ref_label := "p10"],
    pick_representative_cities(S, metric, ref = "mean")[, ref_label := "mean"],
    pick_representative_cities(S, metric, ref = 0.9)[, ref_label := "p90"]
  )
  print(reps)

  # 3) Per-city plots for each representative
  for (ct in reps$city_tag) {
    nf <- file.path(NODES_DIR, paste0("nodes_city_", ct, ".csv"))
    ef <- file.path(EDGES_DIR, paste0("edges_city_", ct, ".csv"))

    if (!file.exists(nf) || !file.exists(ef)) {
      message("Missing files for ", ct, " (skipping per-city plots)")
      next
    }

    nodes <- fread(nf)
    edges <- fread(ef)
    gs <- build_graphs(nodes, edges)

    # Degree distribution
    print(plot_degree_distribution(gs$g_multi, gs$g_simple,
                                   paste0(ct, " | ", metric)))

    # Network map (if coordinates exist)
    if (all(c("longitude", "latitude") %in% names(nodes)) && vcount(gs$g_simple) > 0) {
      plot_network_on_map(gs$g_simple, nodes, paste0(ct, " | ", metric))
      plot_network_on_map(gs$g_multi,  nodes, paste0(ct, " | ", metric))
    }
  }
}

dev.off()

# --- 5c. Cute plots: n_nodes distribution + reference city graphs ---
# Publication-quality plots saved to cute_plots_out/

CUTE_OUT <- "cute_plots_out"
dir.create(CUTE_OUT, showWarnings = FALSE, recursive = TRUE)

res <- readRDS("analytics_out/citylines_metrics.rds")
S   <- res$per_city_summary

# Apply same avg-degree filter
unfilter_n_cities <- nrow(S)
deg_col <- if ("avg_degree_simple" %in% names(S)) "avg_degree_simple"
           else if ("avg_degree_multi" %in% names(S)) "avg_degree_multi"
           else stop("No avg_degree column found.")

S <- S[!is.na(get(deg_col)) & get(deg_col) >= 1.0]
filtered_n_cities <- nrow(S)
message("Cities kept after degree filter (", deg_col, " >= 1): ",
        filtered_n_cities, " (", round(filtered_n_cities / unfilter_n_cities * 100), "%)")

if (nrow(S) == 0) stop("No cities left after filtering by average degree >= 1.")
if (!("n_nodes" %in% names(S))) stop("Column `n_nodes` not found in per_city_summary.")

# Compute n_nodes distribution stats
x <- S[!is.na(n_nodes), n_nodes]
if (length(x) == 0) stop("No non-NA n_nodes values after filtering.")

q10 <- as.numeric(quantile(x, probs = 0.10, names = FALSE, na.rm = TRUE))
q90 <- as.numeric(quantile(x, probs = 0.90, names = FALSE, na.rm = TRUE))
mu  <- mean(x, na.rm = TRUE)

# Choose 3 reference cities (closest to q10, mean, q90)
ref10 <- pick_representative_cities(S, "n_nodes", ref = 0.10)
refmu <- pick_representative_cities(S, "n_nodes", ref = "mean")
ref90 <- pick_representative_cities(S, "n_nodes", ref = 0.90)

ref_cities <- rbindlist(list(ref10, refmu, ref90), use.names = TRUE, fill = TRUE)
ref_cities[, ref := factor(ref, levels = c("p10", "mean", "p90"))]

h         <- hist(x, breaks = 35, plot = FALSE)
max_count <- max(h$counts)
y_cap     <- 15

# PDF 1: Across-cities distribution of n_nodes (no title, publication-ready)
p <- ggplot(data.table(n_nodes = x), aes(x = n_nodes)) +
  geom_histogram(bins = 35) +
  geom_vline(xintercept = q10, linetype = "dashed") +
  geom_vline(xintercept = mu,  linetype = "dashed") +
  geom_vline(xintercept = q90, linetype = "dashed") +
  annotate("text", x = q10, y = y_cap - 0.4, label = "q10",  angle = 90, vjust = 1) +
  annotate("text", x = mu,  y = y_cap - 0.4, label = "mean", angle = 90, vjust = 1) +
  annotate("text", x = q90, y = y_cap - 0.4, label = "q90",  angle = 90, vjust = 1) +
  geom_segment(
    data = ref_cities,
    aes(x = metric_value, xend = metric_value, y = 0, yend = y_cap * 0.06),
    inherit.aes = FALSE, color = "red", linewidth = 0.6
  ) +
  geom_text(
    data = ref_cities,
    aes(x = metric_value, y = 0, label = city_name),
    inherit.aes = FALSE, color = "red", angle = 90, vjust = 1.4
  ) +
  labs(title = NULL, x = "n nodes", y = "cities") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  theme_minimal(base_size = 16) +
  theme(plot.margin = margin(10, 10, 35, 10)) +
  coord_cartesian(ylim = c(0, y_cap), clip = "off")

pdf(file.path(CUTE_OUT, "n_nodes_distribution.pdf"), width = 8.5, height = 4.5, onefile = FALSE)
print(p)
dev.off()

# PDF 2: 3 representative cities (simple graph), horizontally stacked
pdf(file.path(CUTE_OUT, "reference_cities_simple_graphs.pdf"),
    width = 13.5, height = 4.5, onefile = TRUE)
op <- par(mfrow = c(1, 3), mar = c(3.5, 0.5, 0.5, 0.5))

for (i in seq_len(nrow(ref_cities))) {
  ct <- ref_cities$city_tag[i]
  cn <- ref_cities$city_name[i]

  nf <- file.path(NODES_DIR, paste0("nodes_city_", ct, ".csv"))
  ef <- file.path(EDGES_DIR, paste0("edges_city_", ct, ".csv"))

  if (!file.exists(nf) || !file.exists(ef)) {
    plot.new()
    title(main = paste0(cn, "\n(missing files)"))
    next
  }

  nodes <- fread(nf)
  edges <- fread(ef)
  gs    <- build_graphs(nodes, edges)

  if (vcount(gs$g_simple) == 0) {
    plot.new()
    title(main = paste0(cn, "\n(empty graph)"))
    next
  }

  plot_network_on_map(gs$g_simple, nodes, city_tag = "")
  mtext(cn, side = 1, line = 1.2, cex = 2, font = 2, col = "red")
}

par(op)
dev.off()

message("Cute plots written to: ", CUTE_OUT)


# ==============================================================================
# PART 6 - Global Plots for LaTeX Appendix
#
# For each numeric metric (filtered cities, mean degree >= 1):
#   - Compute mean, q10, q90
#   - Pick representative cities closest to q10/mean/q90
#   - Save ONE PDF page per metric:
#       [ histogram + dashed q10/mean/q90 + red rep-city ticks/names ]
#       [ 3 geo network plots (one per representative city) ]
# Output folder: global_plots/
# ==============================================================================
message("\n", strrep("=", 70))
message("PART 6: Global Plots for LaTeX Appendix")
message(strrep("=", 70))

GLOBAL_OUT <- "global_plots"
dir.create(GLOBAL_OUT, showWarnings = FALSE, recursive = TRUE)

# --- Helper: safe filename from metric name ---
.safe_fname <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# Histogram with q10/mean/q90 dashed lines and red representative-city ticks
plot_metric_hist_with_refs <- function(S, metric, reps_dt, q10, mu, q90, bins = 35) {
  x <- S[[metric]]
  x <- x[!is.na(x)]

  h     <- hist(x, breaks = bins, plot = FALSE)
  y_max <- max(h$counts)
  y_top <- y_max * 1.10

  hist(x, breaks = bins, main = "", xlab = metric, ylab = "cities",
       xaxt = "n", ylim = c(0, y_top))
  axis(1)

  # Dashed reference lines
  abline(v = q10, lty = 2)
  abline(v = mu,  lty = 2)
  abline(v = q90, lty = 2)

  dx <- 0.01 * diff(par("usr")[1:2])
  text(q10 - dx, y_top * 0.98, "q10",  srt = 90, adj = c(1, 0.5))
  text(mu  - dx, y_top * 0.98, "mean", srt = 90, adj = c(1, 0.5))
  text(q90 - dx, y_top * 0.98, "q90",  srt = 90, adj = c(1, 0.5))

  # Red ticks + names for representative cities
  if (nrow(reps_dt) > 0) {
    at_vals <- reps_dt$metric_value
    axis(1, at = at_vals, labels = FALSE, col.ticks = "red", tck = -0.03)
    for (i in seq_len(nrow(reps_dt))) {
      text(
        x = reps_dt$metric_value[i],
        y = par("usr")[3],
        labels = reps_dt$city_name[i],
        col = "red", srt = -30, adj = c(0, 1),
        xpd = NA, cex = 1.5, font = 2
      )
    }
  }
}

# Plot a single city's simple graph in one panel (geo layout, community colors)
plot_city_simple_graph_panel <- function(city_tag, city_name) {
  nf <- file.path(NODES_DIR, paste0("nodes_city_", city_tag, ".csv"))
  ef <- file.path(EDGES_DIR, paste0("edges_city_", city_tag, ".csv"))

  if (!file.exists(nf) || !file.exists(ef)) {
    plot.new()
    mtext(paste0(city_name, "\n(missing files)"), side = 3, line = -1, cex = 1.0, font = 2)
    return(invisible(NULL))
  }

  nodes <- fread(nf)
  edges <- fread(ef)
  gs    <- build_graphs(nodes, edges)

  if (vcount(gs$g_simple) == 0) {
    plot.new()
    mtext(paste0(city_name, "\n(empty graph)"), side = 3, line = -1, cex = 1.0, font = 2)
    return(invisible(NULL))
  }

  plot_network_on_map(gs$g_simple, nodes, city_tag = "")
  title(main = "")
  mtext(city_name, side = 1, line = 1.2, cex = 1.1, font = 2)
}

# --- Load summaries and filter ---
res <- readRDS("analytics_out/citylines_metrics.rds")
S   <- res$per_city_summary

unfilter_n_cities <- nrow(S)

deg_col <- if ("avg_degree_simple" %in% names(S)) "avg_degree_simple"
           else if ("avg_degree_multi" %in% names(S)) "avg_degree_multi"
           else stop("No avg_degree column found.")

S <- S[!is.na(get(deg_col)) & get(deg_col) >= 1.0]
filtered_n_cities <- nrow(S)
message("Cities kept after degree filter (", deg_col, " >= 1): ",
        filtered_n_cities, " (", round(filtered_n_cities / unfilter_n_cities * 100), "%)")

if (nrow(S) == 0) stop("No cities left after filtering by average degree >= 1.")

# All numeric metric columns (exclude identifiers)
id_cols     <- c("city_tag", "city_index", "city_name")
metric_cols <- setdiff(names(S), id_cols)
metric_cols <- metric_cols[sapply(S[, ..metric_cols], is.numeric)]

metric_cols <- metric_cols[
  sapply(metric_cols, function(m) {
    x <- S[[m]]
    x <- x[!is.na(x)]
    length(x) >= 5 && length(unique(x)) >= 2
  })
]

message("Will generate global pages for metrics: ", paste(metric_cols, collapse = ", "))

# --- Loop: one PDF per metric ---
for (metric in metric_cols) {
  x <- S[[metric]]
  x <- x[!is.na(x)]
  if (length(x) < 5) next

  q10 <- as.numeric(quantile(x, 0.10, na.rm = TRUE, names = FALSE))
  q90 <- as.numeric(quantile(x, 0.90, na.rm = TRUE, names = FALSE))
  mu  <- mean(x, na.rm = TRUE)

  # 3 representative cities (closest to q10, mean, q90)
  rep10 <- pick_representative_cities(S, metric, ref = 0.10)
  repmu <- pick_representative_cities(S, metric, ref = "mean")
  rep90 <- pick_representative_cities(S, metric, ref = 0.90)

  reps <- rbindlist(list(rep10, repmu, rep90), use.names = TRUE, fill = TRUE)
  reps <- reps[!is.na(city_tag) & !is.na(metric_value)]

  out_pdf <- file.path(GLOBAL_OUT, paste0(.safe_fname(metric), ".pdf"))
  pdf(out_pdf, width = 18, height = 4.8, onefile = TRUE)

  # Layout: 1 row, 4 cols; histogram wider
  layout(matrix(1:4, nrow = 1), widths = c(1.65, 1, 1, 1))

  # Panel 1: distribution histogram
  par(mar = c(7, 4.2, 1.2, 1.0))
  plot_metric_hist_with_refs(S, metric, reps, q10, mu, q90, bins = 35)

  # Panels 2-4: representative city graph panels
  par(mar = c(3.0, 0.6, 0.6, 0.6))
  if (nrow(reps) >= 1) plot_city_simple_graph_panel(reps$city_tag[1], reps$city_name[1]) else plot.new()
  if (nrow(reps) >= 2) plot_city_simple_graph_panel(reps$city_tag[2], reps$city_name[2]) else plot.new()
  if (nrow(reps) >= 3) plot_city_simple_graph_panel(reps$city_tag[3], reps$city_name[3]) else plot.new()

  dev.off()
  message("Wrote: ", out_pdf)
}

message("\n", strrep("=", 70))
message("Pipeline complete.")
message(strrep("=", 70))
