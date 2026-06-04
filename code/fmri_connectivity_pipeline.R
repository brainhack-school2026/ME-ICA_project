# ============================================================
# STEP 1 — INSTALL & LOAD PACKAGES
# ============================================================
install.packages("pheatmap")
install.packages("signal")

library(RNifti)
library(jsonlite)
library(pheatmap)
library(signal)
library(tidyverse)


# ============================================================
# STEP 2 — SET PATHS
# ============================================================
data_dir <- "../data/ds007661/sub-01/func"

bold_path <- file.path(data_dir,
                       "sub-01_task-handgrasp_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
mask_path <- file.path(data_dir,
                       "sub-01_task-handgrasp_run-1_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz")
conf_path <- file.path(data_dir,
                       "sub-01_task-handgrasp_run-1_desc-confounds_timeseries.tsv")


# ============================================================
# STEP 3 — DOWNLOAD SCHAEFER 200 ATLAS (2 mm MNI space)
# ============================================================
atlas_path <- "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"

if (!file.exists(atlas_path)) {
  download.file(
    url = paste0(
      "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/",
      "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/",
      "Parcellations/MNI/",
      "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
    ),
    destfile = atlas_path,
    mode = "wb"
  )
  message("Schaefer 200 atlas downloaded.")
}


# ============================================================
# STEP 4 — RESAMPLE ATLAS TO MATCH BOLD DIMENSIONS
# Your BOLD is 62x77x82, atlas is 2mm MNI — they won't match
# so we use nearest-neighbour resampling (critical for label images)
# ============================================================
resample_atlas_to_bold <- function(atlas_path, bold_path) {
  atlas <- readNifti(atlas_path)
  bold  <- readNifti(bold_path)

  atlas_arr <- as.array(atlas)
  bold_dims <- dim(bold)[1:3]

  cat("Atlas dims:", dim(atlas_arr)[1:3], "\n")
  cat("Target BOLD dims:", bold_dims, "\n")

  scale <- bold_dims / dim(atlas_arr)[1:3]
  out   <- array(0L, dim = bold_dims)

  for (x in seq_len(bold_dims[1])) {
    ax <- min(round(x / scale[1]), dim(atlas_arr)[1])
    for (y in seq_len(bold_dims[2])) {
      ay <- min(round(y / scale[2]), dim(atlas_arr)[2])
      for (z in seq_len(bold_dims[3])) {
        az           <- min(round(z / scale[3]), dim(atlas_arr)[3])
        out[x, y, z] <- atlas_arr[ax, ay, az]
      }
    }
  }

  result   <- asNifti(out, template = bold)
  out_path <- "Schaefer200_resampled.nii.gz"
  writeNifti(result, out_path)
  cat("Saved resampled atlas to:", out_path, "\n")
  return(out_path)
}

atlas_resampled_path <- resample_atlas_to_bold(atlas_path, bold_path)

# Verify dims now match
cat("Resampled atlas dims:", dim(as.array(readNifti(atlas_resampled_path)))[1:3], "\n")
# Should print: 62 77 82


# ============================================================
# STEP 5 — GET TR FROM JSON SIDECAR
# ============================================================
json_path <- file.path(data_dir,
                       "sub-01_task-handgrasp_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.json")

TR <- fromJSON(json_path)$RepetitionTime
cat("TR is:", TR, "seconds\n")


# ============================================================
# STEP 6 — EXTRACT ROI TIME SERIES
# Averages the BOLD signal within each of the 200 brain parcels
# ============================================================
extract_roi_timeseries <- function(bold_path, atlas_path, mask_path) {

  message("Loading BOLD (this may take a moment)...")
  bold_arr  <- as.array(readNifti(bold_path))
  atlas_vec <- as.vector(as.array(readNifti(atlas_path)))
  mask_vec  <- as.vector(as.array(readNifti(mask_path))) > 0

  dims     <- dim(bold_arr)
  n_t      <- dims[4]
  bold_mat <- matrix(bold_arr, nrow = prod(dims[1:3]), ncol = n_t)

  bold_mat[!mask_vec, ] <- NA

  roi_ids <- sort(unique(atlas_vec[atlas_vec > 0]))
  roi_ts  <- matrix(NA, nrow = n_t, ncol = length(roi_ids),
                    dimnames = list(NULL, paste0("ROI_", roi_ids)))

  for (i in seq_along(roi_ids)) {
    idx <- which(atlas_vec == roi_ids[i] & mask_vec)
    if (length(idx) == 0) next
    roi_ts[, i] <- colMeans(bold_mat[idx, , drop = FALSE], na.rm = TRUE)
  }

  message("Extracted ", ncol(roi_ts), " ROI time series, ", n_t, " timepoints each.")
  return(roi_ts)
}


# ============================================================
# STEP 7 — REGRESS OUT NOISE (WM + CSF)
# ============================================================
regress_confounds <- function(roi_ts, conf_path) {

  conf <- read.table(conf_path, header = TRUE, sep = "\t", na.strings = "n/a")

  keep_cols <- intersect(c("white_matter", "csf"), names(conf))

  if (length(keep_cols) == 0) {
    message("No WM/CSF confounds found — skipping regression.")
    return(roi_ts)
  }

  conf_mat <- as.matrix(conf[, keep_cols, drop = FALSE])
  conf_mat[is.na(conf_mat)] <- 0

  cleaned <- apply(roi_ts, 2, function(y) {
    y[is.na(y)] <- mean(y, na.rm = TRUE)
    residuals(lm(y ~ conf_mat))
  })

  message("Confound regression done.")
  return(cleaned)
}


# ============================================================
# STEP 8 — BANDPASS FILTER (keep 0.01–0.1 Hz brain rhythms)
# ============================================================
bandpass_filter <- function(ts_mat, TR, low = 0.01, high = 0.1) {
  nyq <- (1 / TR) / 2
  bf  <- butter(3, c(low / nyq, high / nyq), type = "pass")

  filtered <- apply(ts_mat, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    as.numeric(filtfilt(bf, x))
  })

  message("Bandpass filter applied (", low, "–", high, " Hz).")
  return(filtered)
}


# ============================================================
# STEP 9 — COMPUTE CONNECTIVITY MATRIX
# Pearson correlation + Fisher's z-transform
# ============================================================
compute_connectivity <- function(roi_ts_clean) {
  r_mat <- cor(roi_ts_clean, use = "pairwise.complete.obs")
  z_mat <- atanh(r_mat)
  diag(z_mat) <- 0
  return(z_mat)
}


# ============================================================
# STEP 10 — RUN THE FULL PIPELINE FOR SUB-01, RUN-1
# Note: uses atlas_resampled_path from Step 4
# ============================================================
roi_ts    <- extract_roi_timeseries(bold_path, atlas_resampled_path, mask_path)
roi_clean <- regress_confounds(roi_ts, conf_path)
roi_filt  <- bandpass_filter(roi_clean, TR)
conn_mat  <- compute_connectivity(roi_filt)

cat("Connectivity matrix dimensions:", dim(conn_mat), "\n")  # should be 200 x 200


# ============================================================
# STEP 11 — PLOT THE MATRIX
# ============================================================
pheatmap(
  mat          = conn_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-1, 1, length.out = 101),
  main         = "sub-01 run-1 Connectivity Matrix (Schaefer 200 ROIs)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize     = 8
)


# ============================================================
# STEP 12 — LOOP ALL 8 SUBJECTS × 2 RUNS, THEN AVERAGE RUNS
# ============================================================
subjects <- paste0("sub-0", 1:8)
runs     <- c("run-1", "run-2")
base_dir <- "/path/to/ds007661/derivatives"   # <-- edit this

all_conn      <- list()
all_conn_runs <- list()

for (sub in subjects) {
  func_dir <- file.path(base_dir, sub, "func")
  run_mats <- list()

  for (run in runs) {
    bold <- file.path(func_dir,
                      sprintf("%s_task-handgrasp_%s_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz",
                              sub, run))
    mask <- file.path(func_dir,
                      sprintf("%s_task-handgrasp_%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz",
                              sub, run))
    conf <- file.path(func_dir,
                      sprintf("%s_task-handgrasp_%s_desc-confounds_timeseries.tsv",
                              sub, run))

    if (!file.exists(bold)) { message("Skipping (not found): ", bold); next }

    ts      <- extract_roi_timeseries(bold, atlas_resampled_path, mask)
    cleaned <- regress_confounds(ts, conf)
    filt    <- bandpass_filter(cleaned, TR)
    mat     <- compute_connectivity(filt)

    run_mats[[run]] <- mat
    message("✓ Done: ", sub, " ", run)
  }

  all_conn_runs[[sub]] <- run_mats

  if (length(run_mats) == 2) {
    all_conn[[sub]] <- (run_mats[["run-1"]] + run_mats[["run-2"]]) / 2
  } else if (length(run_mats) == 1) {
    all_conn[[sub]] <- run_mats[[1]]
  }
}

# Step before 13 you need to save the matrix for each subject per run to run the consistency
# For change subject change the number and the number should 1 for run1 an 2 for run2 --> conn_mat holds the actual \(200 x 200\) functional connectivity matrix you just generated
all_conn_runs[["sub-02"]][[1]] <- conn_mat

# Run this also to see if its organized properly to do the consistency test (r value)
str(all_conn_runs[["sub-02"]])


# ============================================================
# STEP 13 — CHECK RUN-TO-RUN CONSISTENCY
# ============================================================
cat("\n--- Run-to-run consistency (Pearson r of upper triangle) ---\n")

for (sub in subjects) {
  mats <- all_conn_runs[[sub]]
  if (length(mats) < 2) next
  idx <- upper.tri(mats[[1]])
  r   <- cor(mats[[1]][idx], mats[[2]][idx])
  cat(sub, ": r =", round(r, 3), "\n")
}
# Values > 0.7 suggest the matrix is stable enough for clinical use



# ============================================================
# STEP 12 — MANUALLY STORE EACH SUBJECT + RUN
# After running Step 10 for each subject/run, store conn_mat here.
# Format: all_conn_runs[["sub-XX"]][["run-Y"]] <- conn_mat
#
# Initialize the storage list once at the start (run this only once)
# ============================================================
all_conn_runs <- list() #did this already DO NOT TOUCH AGAIN ONLY DO IT ONCE!!!!!!!!!!!!!!

saveRDS(all_conn_runs, "all_conn_runs.rds")

# --- sub-01 ---
# Run Step 10 with sub-01 run-1 paths, then store:
all_conn_runs[["sub-01"]][["run-1"]] <- conn_mat

# Run Step 10 with sub-01 run-2 paths, then store:
all_conn_runs[["sub-01"]][["run-2"]] <- conn_mat

# --- sub-02 ---
all_conn_runs[["sub-02"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-02"]][["run-2"]] <- conn_mat

# --- sub-03 ---
all_conn_runs[["sub-03"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-03"]][["run-2"]] <- conn_mat

# --- sub-04 ---
all_conn_runs[["sub-04"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-04"]][["run-2"]] <- conn_mat

# --- sub-05 ---
all_conn_runs[["sub-05"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-05"]][["run-2"]] <- conn_mat

# --- sub-06 ---
all_conn_runs[["sub-06"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-06"]][["run-2"]] <- conn_mat

# --- sub-07 ---
all_conn_runs[["sub-07"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-07"]][["run-2"]] <- conn_mat

# --- sub-08 ---
all_conn_runs[["sub-08"]][["run-1"]] <- conn_mat
all_conn_runs[["sub-08"]][["run-2"]] <- conn_mat


# ============================================================
# OPTIONAL — CHECK A SUBJECT IS STORED CORRECTLY
# Run this after storing to make sure both runs are there
# ============================================================
str(all_conn_runs)  # double check everything is there

# Should show: List of 2 (or 1 if only one run stored yet)
#  $ run-1: num [1:200, 1:200] ...
#  $ run-2: num [1:200, 1:200] ...


# ============================================================
# STEP 13 — RUN-TO-RUN CONSISTENCY (Pearson r of upper triangle)
# Only runs for subjects that have BOTH runs stored
# ============================================================
subjects <- paste0("sub-0", 1:8)

cat("\n--- Run-to-run consistency (Pearson r of upper triangle) ---\n")

for (sub in subjects) {
  mats <- all_conn_runs[[sub]]
  if (is.null(mats[["run-1"]]) || is.null(mats[["run-2"]])) {
    cat(sub, ": missing a run — skipping\n")
    next
  }
  idx <- upper.tri(mats[["run-1"]])
  r   <- cor(mats[["run-1"]][idx], mats[["run-2"]][idx])
  cat(sub, ": r =", round(r, 3), "\n")
}
# Values > 0.7 suggest the matrix is stable enough for clinical use


# ============================================================
# STEP 14 — AVERAGE RUNS PER SUBJECT (once both runs are stored)
# ============================================================
all_conn <- list()

for (sub in subjects) {
  mats <- all_conn_runs[[sub]]
  if (is.null(mats)) next
  
  if (!is.null(mats[["run-1"]]) && !is.null(mats[["run-2"]])) {
    all_conn[[sub]] <- (mats[["run-1"]] + mats[["run-2"]]) / 2
    cat(sub, ": averaged across runs\n")
  } else if (!is.null(mats[["run-1"]])) {
    all_conn[[sub]] <- mats[["run-1"]]
    cat(sub, ": only run-1 available\n")
  } else if (!is.null(mats[["run-2"]])) {
    all_conn[[sub]] <- mats[["run-2"]]
    cat(sub, ": only run-2 available\n")
  }
}

saveRDS(all_conn, "all_conn_averaged.rds")

# Loading all the Data from the Matrix to run Steps 13 & 14
library(RNifti)
all_conn_runs <- readRDS("all_conn_runs.rds")






install.packages("ggplot2")
library(ggplot2)

library(pheatmap)
library(ggplot2)

# ============================================================
# PRECENTRAL GYRUS ROI ANALYSIS
# In Schaefer 200 parcels, Somatomotor network ROIs (roughly
# ROI 30-67) cover precentral gyrus / primary motor cortex.
# We target the core precentral gyrus ROIs specifically.
# ============================================================

# These are the Somatomotor network ROI indices in Schaefer 200
# that correspond to precentral gyrus / M1
precentral_rois <- 30:67   # Somatomotor network parcels

roi_names <- paste0("ROI_", precentral_rois)

# ============================================================
# STEP 1 — EXTRACT MEAN CONNECTIVITY OF PRECENTRAL GYRUS
# FOR EACH SUBJECT (how strongly is M1 connected to the rest
# of the brain on average?)
# ============================================================
subjects <- paste0("sub-0", 1:8)

precentral_conn <- sapply(subjects, function(sub) {
  mat  <- all_conn[[sub]]                         # 200 x 200 matrix
  idx  <- which(colnames(mat) %in% roi_names)     # find precentral ROIs
  # Mean connectivity of precentral ROIs to ALL other ROIs
  mean(mat[idx, -idx], na.rm = TRUE)
})

cat("\n--- Mean Precentral Gyrus Connectivity per Subject ---\n")
print(round(precentral_conn, 3))


# ============================================================
# STEP 2 — BAR PLOT: PRECENTRAL CONNECTIVITY PER SUBJECT
# ============================================================
df_bar <- data.frame(
  Subject     = names(precentral_conn),
  Connectivity = as.numeric(precentral_conn)
)

ggplot(df_bar, aes(x = Subject, y = Connectivity, fill = Connectivity)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_hline(yintercept = mean(precentral_conn), linetype = "dashed",
             color = "black", linewidth = 0.8) +
  annotate("text", x = 0.6, y = mean(precentral_conn) + 0.005,
           label = "Group mean", hjust = 0, size = 3.5) +
  labs(
    title    = "Mean Precentral Gyrus (M1) Connectivity per Subject",
    subtitle = "Higher values = stronger motor cortex integration\nRelevant for Parkinson's & Stroke detection via ME-ICA",
    x        = "Subject",
    y        = "Mean Fisher z (connectivity to rest of brain)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")


# ============================================================
# STEP 3 — HEATMAP: PRECENTRAL GYRUS ROI CONNECTIVITY
# ACROSS THE GROUP (which M1 parcels are most connected?)
# ============================================================
group_avg <- Reduce("+", all_conn) / length(all_conn)

idx <- which(colnames(group_avg) %in% roi_names)

precentral_mat <- group_avg[idx, idx]

pheatmap(
  mat           = precentral_mat,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  color         = colorRampPalette(c("blue", "white", "red"))(100),
  breaks        = seq(-1, 1, length.out = 101),
  main          = "Precentral Gyrus (M1) Inter-ROI Connectivity\n(Group Average, Schaefer 200)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize      = 9
)


# ============================================================
# STEP 4 — FIND TOP 10 STRONGEST CONNECTIONS TO M1
# (which brain regions are most connected to motor cortex?)
# ============================================================
non_idx <- which(!colnames(group_avg) %in% roi_names)

# Average connectivity from each non-M1 ROI to all M1 ROIs
conn_to_m1 <- sapply(non_idx, function(j) {
  mean(group_avg[idx, j], na.rm = TRUE)
})

names(conn_to_m1) <- colnames(group_avg)[non_idx]
top10 <- sort(conn_to_m1, decreasing = TRUE)[1:10]

cat("\n--- Top 10 ROIs Most Connected to Precentral Gyrus ---\n")
print(round(top10, 3))

# Plot it
df_top <- data.frame(
  ROI          = factor(names(top10), levels = rev(names(top10))),
  Connectivity = as.numeric(top10)
)

ggplot(df_top, aes(x = ROI, y = Connectivity, fill = Connectivity)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_flip() +
  labs(
    title    = "Top 10 ROIs Most Connected to Precentral Gyrus (M1)",
    subtitle = "Group average | Handgrasp task fMRI | Schaefer 200",
    x        = "ROI",
    y        = "Mean Fisher z to M1"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")


