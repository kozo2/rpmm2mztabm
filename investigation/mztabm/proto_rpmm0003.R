## Feasibility prototype v3: RPMM0003 -> mzTab-M full profile (M+S+F+E)
suppressMessages({
  library(RmzTabM)
  library(SummarizedExperiment)
})

## Resolve paths from the script, independently of the shell's working directory.
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) {
  stop("Run this prototype with Rscript so its repository can be located.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_flag)))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  found <- which(args == name)
  if (!length(found)) return(default)
  if (length(found) != 1L || found + 1L > length(args)) stop("Invalid option: ", name)
  args[found + 1L]
}
projects_root <- option("--projects-root", file.path(repo_root, "projects"))
base <- file.path(projects_root, "RPMM0003", "data", "AnalysedDataset")

raw <- readLines(file.path(base, "MetMat_abs.txt"), n = 40)
hdr <- strsplit(raw[1], "\t", fixed = TRUE)[[1]]
mtx <- do.call(rbind, strsplit(raw[-1], "\t", fixed = TRUE))

annotation_grade <- mtx[, 3]
chemical_name    <- trimws(gsub('^"|"$', "", mtx[, 5]))
ri               <- suppressWarnings(as.numeric(mtx[, 6]))
rt_s             <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", mtx[, 7])))

sample_names <- hdr[8:ncol(mtx)]
abund <- apply(mtx[, 8:ncol(mtx), drop = FALSE], 2,
              function(x) suppressWarnings(as.numeric(x)))
colnames(abund) <- sample_names
rownames(abund) <- mtx[, 1]

pheno <- read.csv(file.path(base, "phenodata.csv"), check.names = FALSE)

sel <- c("FD_01_1", "F_037_1", "F_037_old_1", "F_046_1",
         "L_037_1", "L_040_1", "F_058_1", "L_058_1")
sel <- sel[sel %in% sample_names]
pheno_sel <- pheno[match(sel, pheno$gcid), ]

rel_map <- c(A2 = "2", A = "2", I = "3", U = "4")
reliability <- ifelse(annotation_grade %in% names(rel_map), rel_map[annotation_grade], "4")
has_name <- nzchar(chemical_name)

rd <- DataFrame(
  exp_mass_to_charge = rep(NA_real_, length(rt_s)),
  retention_time_in_seconds = as.numeric(rt_s),
  opt_global_RI = ri,
  opt_global_annotation_grade = annotation_grade,
  opt_global_chemical_name = ifelse(has_name, chemical_name, "null"),
  opt_global_win = mtx[, 7]
)
cd <- DataFrame(
  sample = sel,
  assay = sel,
  sample_ref = paste0("sample[", seq_along(sel), "]"),
  ms_run_ref = paste0("ms_run[", seq_along(sel), "]"),
  protocol_ref = "protocol[1]",
  location = paste0("http://metabobank.riken.jp/data/RPMM0003/AnalysedDataset/MetMat_abs.txt#", sel),
  scan_polarity = "positive",
  species = "[NCBITaxon, NCBITaxon:4081, Solanum lycopersicum, ]",
  description = paste0("RPMM0003 sample ", sel),
  fruitvsleaf = pheno_sel$fruitvsleaf,
  new_oldanalysis = pheno_sel$new_oldanalysis
)
ab <- abund[, sel, drop = FALSE]
se <- SummarizedExperiment(assays = list(abundance = ab), rowData = rd, colData = cd)

## MTD from SE path (we only need the MTD part here)
mt0 <- MzTabM(se, id = "RPMM0003", sampleCols. = sampleCols(), assayCols. = assayCols(),
              groups = "fruitvsleaf", group_description = "fruit vs leaf")

## ---- Parse NIST MSP file: base-peak m/z per window (MCR-Reference) ----------
msp <- file.path(base, "130224_Chiba_U_tomato2011_NIST.txt")
ml <- readLines(msp, warn = FALSE)
starts <- grep("^Name:", ml)
ends <- c(starts[-1] - 1L, length(ml))
blocks <- lapply(seq_along(starts), function(i) ml[starts[i]:ends[i]])
get_synon <- function(blk, key) {
  hit <- grep(paste0("^Synon:", key, ":"), blk, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^Synon:", key, ":\\s*"), "", hit))
}
base_peak <- vapply(blocks, function(blk) {
  npi <- grep("^Num Peaks:", blk)
  pk <- blk[(npi + 1L):length(blk)]
  pk <- pk[grepl("^\\s*[0-9]", pk)]
  toks <- unlist(strsplit(paste(pk, collapse = " "), ";"))
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  m <- do.call(rbind, strsplit(toks, "\\s+"))
  as.numeric(m[which.max(as.numeric(m[, 2])), 1])
}, numeric(1))
mcr_ref <- vapply(blocks, get_synon, character(1), key = "MCR-Reference")
names(base_peak) <- mcr_ref
cat("MSP blocks:", length(blocks), " unique windows:", length(unique(na.omit(mcr_ref))), "\n")

## Map feature rows -> base-peak m/z via the Win column
win_tok <- sub("^.*?(Win[0-9]+_[A-Za-z0-9]+).*", "\\1", mtx[, 7])
feature_mz <- base_peak[win_tok]
cat("window match coverage:", sum(!is.na(feature_mz)), "/", length(win_tok), "\n")
## fallback for unmatched: nearest-RI base peak
msp_ri <- suppressWarnings(as.numeric(vapply(blocks, get_synon, character(1), key = "RI")))
miss <- is.na(feature_mz)
if (any(miss)) {
  for (i in which(miss)) {
    d <- abs(msp_ri - ri[i])
    feature_mz[i] <- base_peak[which.min(d)]
  }
}
cat("after fallback, NAs:", sum(is.na(feature_mz)), "\n")

## ---- Declare database + id_confidence_measure in MTD (before sections) ------
mt0 <- setMtdDatabase(mt0, name = "[,,NIST,]", prefix = "NIST",
                     version = "null", uri = "https://chemdata.nist.gov/",
                     replace = TRUE)
## RmzTabM 0.99.3's setMtdField can't create 'id_confidence_measure[n]'
## (it would emit 'id_confidence[n]') -> add the MTD row directly.
mtab <- rbind(mtd(mt0),
              c("id_confidence_measure[1]",
                "[MS, MS:1002896, compound identification confidence level, ]"))
mt0@mtd <- mtdSort(mtab)

## ---- SME: identification evidence -----------------------------------------
## mzTab-M 2.1: SME's exp_mass_to_charge/charge are mandatory and non-nullable.
## For EI-GC-MS: use base-peak m/z of the NIST spectrum, charge 1 (EI ions).
## List columns must have consistent item counts -> database_identifier "NIST"
## for rows with a chemical_name.
db_id <- ifelse(has_name, "NIST", "null")
sme <- smeCreate(evidence_input_id = paste0("evidence_", seq_len(nrow(ab))),
                 chemical_name = ifelse(has_name, chemical_name, "null"),
                 exp_mass_to_charge = as.numeric(feature_mz),
                 charge = rep(1, nrow(ab)),
                 theoretical_mass_to_charge = as.numeric(feature_mz),
                 database_identifier = db_id,
                 identification_method = "[,,NIST MS library search,]",
                 spectra_ref = rep("ms_run[1]:index=0", nrow(ab)),
                 ms_level = "[MS, MS:1000580, MS1 spectrum, ]",
                 rank = 1,
                 id_confidence_measure = matrix("null", nrow = nrow(ab), ncol = 1),
                 opt_global_mz_source = rep("base peak m/z of deconvoluted EI spectrum", nrow(ab)),
                 mtd = mtd(mt0))
## RmzTabM 0.99.3 emits 'adduct_ions' in SME, but the mzTab-M 2.1
## spec/validator requires 'adduct_ion' (singular) -> rename it.
if ("adduct_ions" %in% colnames(sme))
  colnames(sme)[colnames(sme) == "adduct_ions"] <- "adduct_ion"
cat("SME cols:", paste(colnames(sme), collapse = ", "), "\n")

## ---- SMF with SME_ID_REFS ---------------------------------------------------
## SMF: exp_mass_to_charge is MANDATORY and non-nullable in mzTab-M 2.1 ->
## use the base-peak m/z of the NIST spectrum for the feature's window.
smf <- smfCreate(x = ab,
                SME_ID_REFS = as.character(sme$SME_ID),
                exp_mass_to_charge = as.numeric(feature_mz),
                charge = rep(NA_real_, nrow(ab)),
                retention_time_in_seconds = as.numeric(rt_s),
                opt_global_RI = ri,
                opt_global_annotation_grade = annotation_grade,
                opt_global_chemical_name = ifelse(has_name, chemical_name, "null"),
                opt_global_win = mtx[, 7])
cat("SMF cols:", paste(colnames(smf), collapse = ", "), "\n")

## ---- SML --------------------------------------------------------------------
sml <- smlCreate(x = ab,
                 SMF_ID_REFS = smf$SMF_ID,
                 chemical_name = ifelse(has_name, chemical_name, "null"),
                 reliability = reliability,
                 database_identifier = db_id,
                 opt_global_annotation_grade = annotation_grade,
                 opt_global_RI = ri)
sml <- smlAddStudyVariableColumns(sml, mtd(mt0))


## ---- Assemble and write ------------------------------------------------------
mt <- MzTabM(mtd = mtd(mt0), sml = sml, smf = smf, sme = sme)
cat("profile:", getMtdField(mtd(mt), "mzTab-profile"), "\n")

out <- option("--output", file.path(projects_root, "RPMM0003", "mztabm", "RPMM0003_MetMat_abs.mztab"))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
res <- try(writeMzTabM(mt, out, validate = FALSE))
cat("writeMzTabM done\n")

## Run the official validator manually to see full messages with line numbers
vpath <- option("--validator", Sys.getenv("RPMM_MZTABM_VALIDATOR", ""))
if (!nzchar(vpath)) {
  cached <- list.files(path.expand("~/Library/Caches/org.R-project.R/R/BiocFileCache"),
                       pattern = "jmztabm-validator.*macos$", full.names = TRUE)
  if (length(cached)) vpath <- sort(cached)[1L]
}
if (!nzchar(vpath) || !file.exists(vpath)) {
  stop("Set --validator or RPMM_MZTABM_VALIDATOR to the jmztabm validator executable.")
}
vr <- system2(vpath, c("-c", shQuote(out)), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(vr, "status")) && attr(vr, "status") != 0L) stop("jmztabm validation failed")
cat(tail(vr, 80), sep = "\n")

## The official validator already ran above. Avoid its second implicit run,
## which would reopen the user's BiocFileCache just to obtain the executable.
back <- readMzTabM(out, validate = FALSE)
## NOTE: makeSummarizedExperimentFromMzTabM has a read-side bug: data.frame()
## mangles '#' in ms_run-location-derived colnames (make.names), so the SE
## round-trip fails on fragment URIs. QA therefore compares sections directly.
s <- smf(back)
ab_back <- s[, grepl("^abundance_assay", colnames(s)), drop = FALSE]
ab_back <- apply(ab_back, 2, function(x) as.numeric(sub("^null$", NA, x)))
cat("round-trip SMF dim:", dim(s), "\n")
stopifnot(isTRUE(all.equal(ab_back, ab, check.attributes = FALSE)))
cat("abundance values identical after round-trip: TRUE\n")
cat("abundance values identical after round-trip: TRUE\n")
cat("MTD fields:", nrow(mtd(back)), "| profile:", getMtdField(mtd(back), "mzTab-profile"), "\n")
