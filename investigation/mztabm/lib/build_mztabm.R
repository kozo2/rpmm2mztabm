## build_mztabm.R - assemble an RmzTabM MzTabM object from a parsed matrix
##
## Ports the verified prototype (proto_rpmm0003.R) to arbitrary matrices and
## centralizes every RmzTabM 0.99.3 workaround catalogued in
## ../mztabm-conversion-plan.md section 5.

suppressMessages({
  library(RmzTabM)
  library(SummarizedExperiment)
})

#' CV parameter helper: "[accession-prefix, accession, name, ]"
cv <- function(prefix, acc, name) paste0("[", prefix, ", ", acc, ", ", name, ", ]")
null_cv <- "[, , null, null]"

#' Reliability mapping from RPMM annotation grade to mzTab-M levels 1-4.
grade_to_reliability <- function(grade) {
  u <- toupper(trim(grade))
  ifelse(grepl("^A", u), "2",
  ifelse(grepl("^I", u), "3",
  ifelse(grepl("^M", u), "2",     # "M2"/"M1" model-predicted annotations -> level 2 with level noted in opt
  ifelse(grepl("^U", u), "4",
  ifelse(nzchar(u), "4", "4")))))
}

#' Build and write one mzTab-M file.
#'
#' @param mx        list from read_matrix()
#' @param pheno     list from read_phenodata() or NULL
#' @param msp       data.frame from parse_msp() or NULL
#' @param species   named character: sample or "" -> species free text / CV string
#' @param cfg       list: project, project_title, species_cv, tissue_cv,
#'                  polarity ("positive"/"negative"), protocol, database_* ,
#'                  description, publication
#' @param named_lut named character: metabolite name -> molecular formula
#' @param out_path  output .mztab file
#' @return list(ok, error, n_features, n_samples, n_dropped_rows, notes)
build_mztabm <- function(mx, pheno, msp, cfg, named_lut = NULL, out_path,
                         validate = TRUE) {
  notes <- character()
  meta_str <- function(x) {
    if (is.null(x) || length(x) == 0L) return("")
    x <- as.character(x)[1L]
    if (is.na(x) || !nzchar(x)) "" else x
  }
  abund <- mx$abund
  n <- nrow(abund)
  msnames <- mx$sample_names

  ## ---- feature m/z (mandatory, non-nullable) ---------------------------
  if (mx$family == "gcms_mcr") {
    fz <- solve_feature_mz(mx$win_token, mx$ri, msp)
    mz <- fz$mz
    mz_src <- fz$source
    ## RT: prefer matrix RT column; fall back to MSP RT by win token
    rt_s <- mx$rt_s
    if (all(is.na(rt_s)) && !is.null(msp)) {
      idx <- msp_win_index(msp)
      rt_s <- vapply(mx$win_token, function(t) {
        if (!is.na(t) && !is.null(idx[[t]])) idx[[t]]$rt_s[1L] else NA_real_
      }, numeric(1))
    }
    db <- ifelse(!is.na(mx$feature_name) & nzchar(mx$feature_name),
                 "NIST", "null")
    db[db == "NIST" & is.na(mz)] <- "null"
  } else {
    rt_s <- mx$rt_s
    ## named metabolite family: m/z from curated formula LUT
    formulas <- vapply(mx$feature_name, function(nm) {
      if (is.null(named_lut)) return("")
      lkp <- named_lut[tolower(trim(nm))]
      if (is.null(lkp) || !length(lkp) || is.na(lkp)) "" else lkp[1L]
    }, character(1))
    theory_neutral <- vapply(formulas, mass_plus_H, numeric(1))
    mz <- theory_neutral                    # [M+H]+
    mz_src <- ifelse(is.na(mz), NA_character_,
      "theoretical [M+H]+ of named compound (isomer/derivative unresolved)")
    db <- ifelse(!is.na(mz), "RPMM", "null")
    mx$feature_name[is.na(mx$feature_name) | !nzchar(mx$feature_name)] <- NA
    notes <- c(notes, paste0("named-metabolite LUT coverage: ",
                             sum(!is.na(mz)), "/", n))
  }

  ## rows without defensible m/z cannot satisfy mandatory columns -> drop
  keep <- !is.na(mz)
  n_dropped <- sum(!keep)
  dropped_ids <- if (n_dropped) mx$feature_id[!keep] else character()
  if (n_dropped)
    notes <- c(notes, paste0("dropped ", n_dropped,
                             " feature(s) with no defensible m/z"))
  if (!sum(keep)) stop("no feature has a defensible m/z: ", mx$path)
  abund <- abund[keep, , drop = FALSE]
  mz <- mz[keep]; mz_src <- mz_src[keep]; db <- db[keep]
  rt_s <- rt_s[keep]
  feature_name <- mx$feature_name[keep]
  feature_id <- mx$feature_id[keep]
  grade <- mx$grade[keep]
  ri <- mx$ri[keep]
  win <- mx$win_token[keep]
  theory_neutral <- if (mx$family == "gcms_mcr") rep(NA_real_, length(mz))
                    else theory_neutral[keep]
  reliability <- grade_to_reliability(grade)
  has_name <- !is.na(feature_name) & nzchar(feature_name)
  chemical_name <- ifelse(has_name, feature_name, "null")
  n <- nrow(abund)

  ## ---- colData: 1 sample == 1 assay == 1 ms_run ------------------------
  has_pheno <- !is.null(pheno) && length(pheno$ids)
  pheno_row <- if (has_pheno) match(msnames, pheno$ids) else rep(NA_integer_, length(msnames))
  pheno_row <- ifelse(is.na(pheno_row), NA_integer_, pheno_row)
  pheno_get <- function(col, default = NA_character_) {
    out <- rep(default, length(msnames))
    if (has_pheno && col %in% colnames(pheno$df)) {
      v <- pheno$df[[col]][pheno_row]
      v[is.na(v)] <- default
      out <- trim(v)
      out[out == ""] <- default
    }
    out
  }

  ## choose study-variable groups: compact columns preferred, species excluded
  groups <- character()
  if (has_pheno) {
    cand <- pheno$groups
    cand <- cand[!grepl(paste(c("species", "sepices", "speacies", "binomial",
                                "speciesname", "plantspecies", "weight", "runorder",
                                "matlab", "order"), collapse = "|"),
                        tolower(cand))]
    prio <- grepl("class|group|genotype|cultivar|treatment|type|variety|day|time",
                  tolower(cand))
    cand <- c(cand[prio], cand[!prio])
    groups <- utils::head(cand, 2L)
  }

  location_base <- paste0("http://metabobank.riken.jp/data/",
                          sub("investigation/downloads/data/", "",
                              dirname(mx$path)))
  cd <- DataFrame(
    sample = make.unique(msnames),
    assay = make.unique(msnames),
    sample_ref = paste0("sample[", seq_along(msnames), "]"),
    ms_run_ref = paste0("ms_run[", seq_along(msnames), "]"),
    protocol_ref = "protocol[1]",
    location = file.path(location_base, paste0(basename(mx$path), "#", msnames)),
    scan_polarity = cfg$polarity %||% "positive",
    description = paste0(cfg$project, " sample ", msnames))
  if (nzchar(meta_str(cfg$species_cv)))
    cd$species <- cfg$species_cv
  if (nzchar(meta_str(cfg$tissue_cv)))
    cd$tissue <- cfg$tissue_cv
  for (g in groups) {
    v <- pheno_get(g, "NA")
    cd[[g]] <- paste0(g, "=", v)     # unique-readable study variable values
  }

  se <- SummarizedExperiment(assays = list(abundance = abund),
                             rowData = DataFrame(row_id = feature_id),
                             colData = cd)
  mt0 <- MzTabM(se, id = cfg$id, sampleCols. = sampleCols(),
                assayCols. = assayCols(), groups = groups,
                group_description = paste0("grouping variable ", groups))

  ## ---- MTD extras ------------------------------------------------------
  desc <- meta_str(cfg$description)
  if (!nzchar(desc)) desc <- paste("Converted from", basename(mx$path))
  ttl <- meta_str(cfg$title)
  if (!nzchar(ttl)) ttl <- cfg$id
  mt0 <- setMtdField(mt0, "description", desc)
  mt0 <- setMtdField(mt0, "title", ttl)
  if (nzchar(meta_str(cfg$publication)))
    mt0 <- setMtdField(mt0, "publication", cfg$publication)
  if (nzchar(meta_str(cfg$external_study)))
    mt0 <- setMtdField(mt0, "external_study", cfg$external_study)
  mt0 <- setMtdDatabase(mt0, name = "[,,NIST,]", prefix = "NIST",
                        version = "null", uri = "https://chemdata.nist.gov/",
                        replace = FALSE)
  if (any(db == "RPMM"))
    mt0 <- setMtdDatabase(mt0, name = "[,,RPMM,]", prefix = "RPMM", version = "null",
                          uri = "http://metabobank.riken.jp/pmm/", replace = FALSE)

  ## id_confidence_measure (RmzTabM 0.99.3 can't create the [<prefix>, <id>]
  ## indexed field name -> append the MTD row manually)
  mtab <- rbind(mtd(mt0),
                c("id_confidence_measure[1]",
                  "[MS, MS:1002896, compound identification confidence level, ]"))
  mt0@mtd <- mtdSort(mtab)
  ## NOTE: mzTab-M 2.1 protocol-block labels differ from what RmzTabM 0.99.3's
  ## setMtdProtocol() emits ("protocol[n]-name" is rejected by the jmztabm
  ## validator). Until the legal shape is established by probe (see plan §5),
  ## the protocol free text goes into MTD description instead - never invented
  ## schema content.
  protocol_note <- meta_str(cfg$protocol)
  if (nzchar(protocol_note))
    mt0 <- setMtdField(mt0, "description", paste0(desc, " Protocol: ", protocol_note))

  ## ---- SME ------------------------------------------------------------
  evidence_ids <- paste0("evidence_", seq_len(n))
  idmeas <- matrix("null", nrow = n, ncol = 1L)
  sme <- smeCreate(
    evidence_input_id = evidence_ids,
    chemical_name = chemical_name,
    exp_mass_to_charge = as.numeric(mz),
    charge = rep(1, n),
    theoretical_mass_to_charge = as.numeric(mz),
    database_identifier = db,
    identification_method = ifelse(has_name,
      if (mx$family == "gcms_mcr") "[,,NIST MS library search,]"
      else "[,,RPMM reviewed annotation,]", "[, , null, null]"),
    spectra_ref = rep("ms_run[1]:index=0", n),
    ms_level = if (mx$family == "gcms_mcr") "[MS, MS:1000580, MS1 spectrum, ]"
               else "[MS, MS:1002421, unknown, ]",
    rank = rep(1, n),
    id_confidence_measure = idmeas,
    opt_global_mz_source = ifelse(is.na(mz_src), "null", mz_src),
    mtd = mtd(mt0))
  if ("adduct_ions" %in% colnames(sme))
    colnames(sme)[colnames(sme) == "adduct_ions"] <- "adduct_ion"

  ## ---- SMF ------------------------------------------------------------
  smf <- smfCreate(
    x = abund,
    SME_ID_REFS = as.character(sme$SME_ID),
    exp_mass_to_charge = as.numeric(mz),
    charge = rep(NA_real_, n),
    retention_time_in_seconds = as.numeric(rt_s),
    adduct_ion = ifelse(mx$family == "gcms_mcr", "null", "[M+H]+"),
    opt_global_RI = ifelse(is.na(ri), "null", format(ri, trim = TRUE)),
    opt_global_annotation_grade = ifelse(is.na(grade) | !nzchar(grade), "null", grade),
    opt_global_chemical_name = chemical_name,
    opt_global_win = ifelse(is.na(win) | !nzchar(win), "null", win),
    opt_global_feature_id = feature_id)

  ## ---- SML ------------------------------------------------------------
  sml <- smlCreate(
    x = abund,
    SMF_ID_REFS = smf$SMF_ID,
    chemical_name = chemical_name,
    reliability = reliability,
    database_identifier = db,
    theoretical_neutral_mass = as.numeric(theory_neutral),
    ## NOTE: SML best_id_confidence_* must be a CV parameter + Double per the
    ## jmztabm validator; RPMM grades ('IIS', 'A1', ...) are neither. The grade
    ## is preserved verbatim in opt_global_annotation_grade instead.
    opt_global_annotation_grade = ifelse(is.na(grade) | !nzchar(grade), "null", grade),
    opt_global_RI = ifelse(is.na(ri), "null", format(ri, trim = TRUE)),
    opt_global_feature_id = feature_id)
  sml <- smlAddStudyVariableColumns(sml, mtd(mt0))

  ## ---- Assemble / write ----------------------------------------------
  mt <- MzTabM(mtd = mtd(mt0), sml = sml, smf = smf, sme = sme)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeMzTabM(mt, out_path, validate = validate)

  list(ok = TRUE, error = NA_character_,
       n_features = n, n_samples = ncol(abund), n_dropped = n_dropped,
       dropped_ids = dropped_ids, kept_abund = abund,
       notes = notes, groups = groups, has_pheno = has_pheno)
}
