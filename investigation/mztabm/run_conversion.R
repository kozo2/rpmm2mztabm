## run_conversion.R - batch-convert RPMM quantification matrices to mzTab-M
##
## Usage:
##   Rscript investigation/mztabm/run_conversion.R [projects ...]
##   Rscript investigation/mztabm/run_conversion.R --inventory-only
##
## Output:
##   investigation/downloads/mztabm/<RPMM_ID>/<RPMM_ID>_<MatrixName>.mztab
##   investigation/mztabm/out/conversion_qa.tsv      agreed at runtime
##   investigation/mztabm/out/curation_issues.tsv     fail-and-log list
##   investigation/mztabm/out/manifest.tsv            the pairing used
##
## Design: everything that cannot be converted faithfully is logged, never
## fabricated (plan section on failure policy). Each conversion is QA-checked
## by re-reading the written file's SMF with readMzTabM() and comparing
## abundance values to the source matrix.

.lib <- "investigation/mztabm/lib"
source(file.path(.lib, "io.R"))
source(file.path(.lib, "read_matrix.R"))
source(file.path(.lib, "parse_msp.R"))
source(file.path(.lib, "inventory.R"))
source(file.path(.lib, "build_mztabm.R"))

DATA   <- "investigation/downloads/data"
MZTABM <- "investigation/downloads/mztabm"
OUT    <- "investigation/mztabm/out"
LUT    <- "investigation/mztabm/data/named_metabolite_formulas.tsv"

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- project titles from the harvested Metabolonote metadata --------------
projects <- local({
  p <- "investigation/projects.json"
  if (!file.exists(p)) return(NULL)
  j <- jsonlite::fromJSON(p, simplifyVector = TRUE)
  ids <- sub(".*/", "", j$p)
  setNames(j$label, ids)
})

## ---- curated formula LUT for named-metabolite matrices --------------------
named_lut <- if (file.exists(LUT)) {
  lu <- read.delim(LUT, stringsAsFactors = FALSE, quote = "")
  lu <- lu[nzchar(lu$metabolite) & nzchar(lu$formula), , drop = FALSE]
  if (!nrow(lu)) NULL else
    setNames(lu$formula, tolower(trim(lu$metabolite)))
} else {
  NULL
}

## ---- QA: re-read the written file and compare abundances -----------------
## Per plan: compare parsed SMF sections directly (not makeSummarizedExperiment
## -- known RmzTabM read-side bug compacts '#' column names).
qa_roundtrip <- function(out_path, kept_abund) {
  tryCatch({
    back <- suppressMessages(readMzTabM(out_path))
    s <- back@smf
    ab_cols <- grep("^abundance_(assay|sample)\\[", colnames(s))
    if (!length(ab_cols)) ab_cols <- grep("^abundance_", colnames(s))
    if (!length(ab_cols)) return(list(ok = FALSE, note = "no abundance columns written"))
    written <- apply(s[, ab_cols, drop = FALSE], 2, function(x) suppressWarnings(as.numeric(as.character(x))))
    if (!identical(dim(written), dim(kept_abund)))
      return(list(ok = FALSE, note = sprintf(
        "dim mismatch: written %dx%d src %dx%d",
        nrow(written), ncol(written), nrow(kept_abund), ncol(kept_abund))))
    eq <- all.equal(unname(written), unname(kept_abund), tolerance = 1e-6,
                    check.attributes = FALSE)
    list(ok = isTRUE(eq),
         note = if (isTRUE(eq)) "abundances identical" else paste(eq, collapse = "; "))
  }, error = function(e) {
    list(ok = FALSE, note = paste("QA read error:", conditionMessage(e)))
  })
}

## ---- CLI ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
inv_only <- "--inventory-only" %in% args
only <- args[!startsWith(args, "--")]

cat("== building inventory ==\n")
inv <- build_inventory(DATA, projects)
write.table(inv$manifest, file.path(OUT, "manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(inv$curation, file.path(OUT, "curation_issues.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("manifest rows:", nrow(inv$manifest),
    " curation issues:", nrow(inv$curation), "\n")
if (inv_only) quit(save = "no")

man <- inv$manifest
if (length(only)) man <- man[man$project %in% only, , drop = FALSE]
if (!nrow(man)) stop("no matrices selected")

## ---- conversion loop ------------------------------------------------------
qa_rows <- list()
for (i in seq_len(nrow(man))) {
  r <- man[i, ]
  stem <- tools::file_path_sans_ext(basename(r$matrix_rel))
  out_path <- file.path(MZTABM, r$project, paste0(r$project, "_", stem, ".mztab"))
  cat(sprintf("[%d/%d] %s -> %s\n", i, nrow(man), r$matrix_rel, out_path))

  row <- data.frame(
    project = r$project, matrix = r$matrix_rel, family = r$family,
    out = out_path, status = "unknown", n_features = NA_integer_,
    n_samples = NA_integer_, n_dropped = NA_integer_, qa = "", note = "",
    stringsAsFactors = FALSE)

  res <- tryCatch({
    mx <- read_matrix(r$matrix_path)
    ph <- if (nzchar(r$phenodata_path))
      tryCatch(read_phenodata(r$phenodata_path, mx$sample_names),
               error = function(e) { row$note <<- paste("phenodata:", conditionMessage(e)); NULL })
    else NULL
    sp <- if (nzchar(r$msp_path)) parse_msp(r$msp_path) else NULL

    cfg <- list(
      id = paste0(r$project, "_", stem),
      project = r$project,
      title = if (nzchar(r$title)) paste0(r$title, " - ", stem) else stem,
      description = paste0("Converted from RPMM dataset file ", r$matrix_rel,
                           " (Metabolonote project ", r$project, ")."),
      publication = character(),
      external_study = paste0("http://metabobank.riken.jp/pmm/", r$project),
      species_cv = r$species_cv,
      tissue_cv = character(),
      polarity = r$polarity,
      protocol = character())

    b <- build_mztabm(mx, ph, sp, cfg, named_lut = named_lut,
                      out_path = out_path, validate = TRUE)
    qa <- qa_roundtrip(out_path, b$kept_abund)
    row$status <- if (qa$ok) "converted" else "converted_qa_fail"
    row$n_features <- b$n_features
    row$n_samples <- b$n_samples
    row$n_dropped <- b$n_dropped
    row$qa <- qa$note
    row$note <- paste(b$notes, collapse = "; ")
    b
  }, error = function(e) {
    row$status <<- "failed"
    row$note <<- conditionMessage(e)
    NULL
  })
  qa_rows[[i]] <- row
}

qa <- do.call(rbind, qa_rows)
write.table(qa, file.path(OUT, "conversion_qa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n== summary ==\n")
print(table(qa$status))
cat("converted:", sum(qa$status == "converted"), "/", nrow(qa), "\n")
for (i in which(qa$status == "failed"))
  cat("FAILED:", qa$matrix[i], "->", substr(qa$note[i], 1, 160), "\n")
