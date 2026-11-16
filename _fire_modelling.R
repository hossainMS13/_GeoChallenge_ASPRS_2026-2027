###############################################################################
# ASPRS GeoChallenge 2026 - burn perimeter models (A: index, B: unsupervised)
# Input: <SITE>_S2_pre.tif and <SITE>_S2_post.tif exported from Earth Engine
###############################################################################
library(terra)
library(sf)

# ================= 1. SITE SETTINGS (change per fire) =================
SITE     <- "pinegulch"
IMG_DIR  <- "C:/Users/mdshossa/OneDrive/OneDrive - Florida International University/GeoChallengeWebsite/GeoChallenge_pinegulch"
REF_FILE <- "C:/Users/mdshossa/OneDrive/OneDrive - Florida International University/GeoChallengeWebsite/bsp_rmwgbd61ef8hugzoosou/mtbs/2020/mtbs_co3933610852620200731_10021021/mtbs_co3933610852620200731_10021021_20190818_20210812_burn_area.shp"
# set REF_FILE <- NULL if no reference
LC_FILE  <- NULL                         # optional ESA WorldCover GeoTIFF

# ----- Path check: try the given path, then a common variant, then ask -----
resolve_path <- function(path, what) {
  if (file.exists(path)) return(path)
  alt <- sub("/OneDrive/OneDrive - ", "/OneDrive - ", path, fixed = TRUE)
  if (file.exists(alt)) return(alt)
  message("Not found: ", path, "\nPlease select ", what, " in the dialog.")
  file.choose()
}

pre_file  <- resolve_path(file.path(IMG_DIR, paste0(SITE, "_S2_pre.tif")),
                          "the PRE-fire image (.tif)")
IMG_DIR   <- dirname(pre_file)
post_file <- resolve_path(file.path(IMG_DIR, paste0(SITE, "_S2_post.tif")),
                          "the POST-fire image (.tif)")
if (!is.null(REF_FILE))
  REF_FILE <- resolve_path(REF_FILE, "the reference perimeter (.shp)")

OUT_DIR <- file.path(IMG_DIR, "outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Pre image :", pre_file, "\nPost image:", post_file,
    "\nReference :", if (is.null(REF_FILE)) "none" else REF_FILE,
    "\nOutputs   :", OUT_DIR, "\n")

# ================= 2. MODEL SETTINGS (same for every fire) =================
MIN_HA     <- 10      # smallest polygon kept, hectares
OPEN_PX    <- 1       # opening radius, pixels (removes speckle)
CLOSE_PX   <- 3       # closing radius, pixels (fills gaps)
THR_MIN    <- 0.10    # guard rails for the automatic threshold
THR_MAX    <- 0.35
N_CLUSTERS <- 5       # Model B: k-means clusters
SEED_PCTL  <- 0.75    # Model B: seeds = top quarter of candidate dNBR
N_SAMPLE   <- 50000   # Model B: pixels used to fit k-means

# ================= 3. LOAD IMAGES =================
bands <- c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B11","B12")

load_s2 <- function(f) {
  r <- rast(f)
  if (nlyr(r) != length(bands))
    stop("Expected 12 bands in ", basename(f), " but found ", nlyr(r))
  names(r) <- bands
  nodata <- (r[["B8A"]] == 0) & (r[["B12"]] == 0)   # exported gaps are 0
  mask(r, nodata, maskvalues = 1) / 10000           # reflectance 0-1
}
pre  <- load_s2(pre_file)
post <- load_s2(post_file)
compareGeom(pre, post)                              # stops if the grids differ
cell_ha <- prod(res(pre)) / 1e4

# ================= 4. INDICES AND EXCLUSION MASK =================
nd <- function(r, a, b) (r[[a]] - r[[b]]) / (r[[a]] + r[[b]])

nbr_pre  <- nd(pre,  "B8A", "B12")
nbr_post <- nd(post, "B8A", "B12")
dnbr     <- nbr_pre - nbr_post
dndvi    <- nd(pre, "B8", "B4")   - nd(post, "B8", "B4")
dnbr2    <- nd(pre, "B11", "B12") - nd(post, "B11", "B12")
dswir    <- post[["B12"]] - pre[["B12"]]
names(dnbr) <- "dNBR"

water    <- (nd(pre, "B3", "B11") > 0) | (nd(post, "B3", "B11") > 0)
burnable <- !water
if (!is.null(LC_FILE)) {
  lc <- project(rast(LC_FILE), dnbr, method = "near")
  burnable <- burnable & !(lc %in% c(40, 50, 80))   # cropland, built-up, water
}
d <- mask(dnbr, burnable, maskvalues = 0)

# ================= 5. HELPER FUNCTIONS =================
# Otsu threshold from a numeric vector
otsu <- function(x, bins = 256) {
  x <- x[is.finite(x)]
  q <- quantile(x, c(0.001, 0.999))
  x <- x[x >= q[1] & x <= q[2]]
  h <- hist(x, breaks = bins, plot = FALSE)
  p <- h$counts / sum(h$counts); m <- h$mids
  w0 <- cumsum(p); mu <- cumsum(p * m); muT <- mu[length(mu)]
  bss <- (muT * w0 - mu)^2 / (w0 * (1 - w0))
  m[which.max(bss)]
}

# Opening then closing on a 0/1 raster
clean <- function(b) {
  b <- ifel(is.na(b), 0, b)
  k <- function(n) matrix(1, 2 * n + 1, 2 * n + 1)
  b <- focal(b, k(OPEN_PX),  "min", na.rm = TRUE)
  b <- focal(b, k(OPEN_PX),  "max", na.rm = TRUE)
  b <- focal(b, k(CLOSE_PX), "max", na.rm = TRUE)
  focal(b, k(CLOSE_PX), "min", na.rm = TRUE)
}

# Polygons from a 0/1 raster; optional seed raster keeps only seeded regions
to_perimeter <- function(bin, label, seed = NULL) {
  p  <- patches(bin, directions = 8, zeroAsNA = TRUE)
  fq <- freq(p)
  keep <- fq$value[fq$count * cell_ha >= MIN_HA]
  if (!is.null(seed)) {
    z <- zonal(seed, p, "max", na.rm = TRUE)
    keep <- intersect(keep, z[[1]][z[[2]] == 1])
  }
  if (length(keep) == 0) stop("No burned region left for ", label)
  p <- ifel(p %in% keep, 1, NA)
  s <- st_as_sf(disagg(as.polygons(p, dissolve = TRUE)))
  s <- st_simplify(s, dTolerance = 1.5 * res(bin)[1])
  s$ha <- as.numeric(st_area(s)) / 1e4
  s
}

# Fill enclosed holes so the perimeter is contiguous
fill_holes <- function(p) {
  s <- st_as_sf(fillHoles(vect(p)))
  s$ha <- as.numeric(st_area(s)) / 1e4
  s
}

save_perimeter <- function(p, label) {
  st_write(p, file.path(OUT_DIR, paste0(SITE, "_perimeter_", label, ".shp")),
           delete_dsn = TRUE, quiet = TRUE)
}

# Overlap with the reference perimeter
validate <- function(perim) {
  pr  <- rasterize(vect(perim), dnbr, field = 1, background = 0)
  ha  <- function(x) global(x, "sum", na.rm = TRUE)[1, 1] * cell_ha
  out <- c(mapped_ha    = ha(pr),
           reference_ha = ha(ref),
           intersection = ha(pr & ref),
           union        = ha(pr | ref),
           commission   = ha(pr & !ref),
           omission     = ha(ref & !pr))
  out["precision"] <- out["intersection"] / out["mapped_ha"]
  out["recall"]    <- out["intersection"] / out["reference_ha"]
  out["IoU"]       <- out["intersection"] / out["union"]
  out
}

# ================= MODEL A: INDEX-BASED =================
thr <- otsu(values(d, na.rm = TRUE))
thr <- min(max(thr, THR_MIN), THR_MAX)
cat("Model A dNBR threshold:", round(thr, 3), "\n")

burned_A    <- clean(d > thr)
perim_A_raw <- to_perimeter(burned_A, "A_index")
perim_A     <- fill_holes(perim_A_raw)
save_perimeter(perim_A, "A_index")

# ================= MODEL B: UNSUPERVISED =================
# B1. Feature stack
feat <- c(dnbr, dndvi, dnbr2, dswir, nbr_post, post[["B8A"]])
names(feat) <- c("dNBR", "dNDVI", "dNBR2", "dSWIR", "NBR_post", "NIR_post")
feat <- mask(feat, burnable, maskvalues = 0)

X  <- values(feat)
ok <- which(complete.cases(X))

# B2. Fit k-means on a standardised random sample
set.seed(42)
idx <- sample(ok, min(N_SAMPLE, length(ok)))
mu  <- colMeans(X[idx, ]); sdv <- apply(X[idx, ], 2, sd)
km  <- kmeans(scale(X[idx, ], mu, sdv), centers = N_CLUSTERS,
              nstart = 10, iter.max = 100)

# B3. Assign every valid pixel to its nearest cluster centre
Xo   <- sweep(sweep(X[ok, ], 2, mu), 2, sdv, "/")
best <- rep(Inf, length(ok)); lab <- integer(length(ok))
for (k in seq_len(N_CLUSTERS)) {
  dk  <- rowSums(sweep(Xo, 2, km$centers[k, ])^2)
  upd <- dk < best
  lab[upd] <- k; best[upd] <- dk[upd]
}
dn_ok <- X[ok, "dNBR"]
rm(Xo, X); gc()

v <- rep(NA_integer_, ncell(dnbr)); v[ok] <- lab
clus <- setValues(rast(dnbr), v)
clus <- mask(focal(clus, 3, "modal", na.rm = TRUE), burnable, maskvalues = 0)

# B4. Automatic labelling: burned = clusters whose mean dNBR exceeds the threshold
cm <- tapply(dn_ok, lab, mean)
cat("Model B mean dNBR per cluster:\n"); print(round(cm, 3))
burned_ids <- as.integer(names(cm)[cm > thr])
if (length(burned_ids) == 0) stop("No cluster passed the burn rule")
cat("Clusters labelled burned:", burned_ids, "\n")

# B5. Seeds: the most severe quarter of candidate pixels
seed_thr <- quantile(dn_ok[lab %in% burned_ids], SEED_PCTL)
cat("Model B seed dNBR threshold:", round(seed_thr, 3), "\n")
cand  <- clus %in% burned_ids
seeds <- (dnbr > seed_thr) & cand

# B6. Grow: keep connected candidate regions that contain a seed
burned_B    <- clean(cand)
perim_B_raw <- to_perimeter(burned_B, "B_unsup", seed = seeds)
perim_B     <- fill_holes(perim_B_raw)
save_perimeter(perim_B, "B_unsup")

# ================= 6. RESULTS =================
if (!is.null(REF_FILE)) {
  ref_v <- project(vect(REF_FILE), crs(dnbr))
  ref   <- rasterize(ref_v, dnbr, field = 1, background = 0)
  
  res_tab <- rbind(A_raw    = validate(perim_A_raw),
                   A_filled = validate(perim_A),
                   B_raw    = validate(perim_B_raw),
                   B_filled = validate(perim_B))
  print(round(res_tab, 3))
  write.csv(res_tab, file.path(OUT_DIR, paste0(SITE, "_validation.csv")))
  
  # Data gaps inside the reference perimeter
  na_in_ref <- global(is.na(dnbr) & ref == 1, "sum")[1, 1] * cell_ha
  cat("Reference area with no data (ha):", round(na_in_ref), "\n")
} else {
  cat("Mapped area (ha)  A:", round(sum(perim_A$ha)),
      " B:", round(sum(perim_B$ha)), "\n")
}

# Map 1: dNBR with perimeters
plot(dnbr, range = c(-0.1, 0.7), main = paste(SITE, "- dNBR with perimeters"))
plot(st_geometry(perim_A), add = TRUE, border = "yellow", lwd = 2)
plot(st_geometry(perim_B), add = TRUE, border = "cyan",   lwd = 2)
if (!is.null(REF_FILE)) plot(ref_v, add = TRUE, border = "red", lwd = 2)

# Map 2: where Model B still misses the reference
if (!is.null(REF_FILE)) {
  pred_B <- rasterize(vect(perim_B), dnbr, field = 1, background = 0)
  plot(ifel(ref == 1 & pred_B == 0, 1, NA), col = "red", legend = FALSE,
       main = paste(SITE, "- Model B omission (red)"))
  plot(ref_v, add = TRUE)
}

# global(is.na(pre[["B12"]]),  "sum")    # gap cells in the pre-fire image
# global(is.na(post[["B12"]]), "sum")    # gap cells in the post-fire image
# library(terra)
# 
# # Identify gap cells (NA = gap)
# pre_gap  <- is.na(pre[["B8A"]])
# post_gap <- is.na(post[["B8A"]])
# 
# # Plot side by side
# par(mfrow = c(1, 2))
# 
# plot(pre_gap,
#      main = "Pre-fire gaps",
#      col = c("lightgreen", "red"),
#      legend = FALSE)
# 
# plot(post_gap,
#      main = "Post-fire gaps",
#      col = c("lightgreen", "red"),
#      legend = FALSE)
# 
# par(mfrow = c(1, 1))

# ---- Raw vs filled: three-panel comparison ----
plot_raw_filled <- function(raw, filled, model_name, zoom = NULL) {
  r_raw <- rasterize(vect(raw),    dnbr, field = 1, background = 0)
  r_fil <- rasterize(vect(filled), dnbr, field = 1, background = 0)
  both  <- ifel(r_fil == 0, NA, r_raw + r_fil)      # 1 = added by filling, 2 = detected
  if (!is.null(zoom)) {
    r_raw <- crop(r_raw, zoom); r_fil <- crop(r_fil, zoom); both <- crop(both, zoom)
  }
  add_ref <- function() if (!is.null(REF_FILE)) plot(ref_v, add = TRUE, border = "black", lwd = 1.5)
  
  op <- par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))
  plot(ifel(r_raw == 1, 1, NA), col = "firebrick", legend = FALSE, axes = FALSE,
       main = paste(model_name, "- raw (detected burn)"));  add_ref()
  plot(ifel(r_fil == 1, 1, NA), col = "firebrick", legend = FALSE, axes = FALSE,
       main = paste(model_name, "- filled (perimeter)"));   add_ref()
  plot(both, col = c("gold", "firebrick"), type = "classes",
       levels = c("added by filling", "detected"), axes = FALSE,
       main = paste(model_name, "- what filling added"));   add_ref()
  par(op)
}

# Whole fire, Model B
plot_raw_filled(perim_B_raw, perim_B, "Model B")

# Same for Model A
plot_raw_filled(perim_A_raw, perim_A, "Model A")