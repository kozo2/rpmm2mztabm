## read_matrix.R - read a quantification matrix and classify its columns
##
## Two dominant RPMM families exist:
##   (1) "MCR feature" matrices (GC-MS, 45+ files): a run of annotation/metadata
##       columns then sample abundance columns. Metadata columns vary in count
##       (7 base, +Marked_peaks/Representative/fData::* extras) and marker names
##       are stable enough to be detectable by normalized name.
##   (2) "named metabolite" matrices (Kusano, carodata, LCMS flavonoid): a single
##       metabolite-name column then sample columns.
##
## Strategy: explicitly gather known metadata columns by normalized name; the
## remaining columns are samples iff their first 60 data rows look numeric
## (>= 50%). This deliberately keeps QC / standard columns as samples and only
## drops text metadata.
##
## Header rows are located by anchor-token search, not by assuming row 1:
## MarkerLynx exports carry a banner line first (RPMM0068/Lipid)

`%||%` <- function(a, b) if (is.null(a) || !length(a) || is.na(a)) b else a

## io.R must already be loaded (see load.R); define trim() fallback:
if (!exists("trim", mode = "function")) {
  trim <- function(x) gsub("^\\s+|\\s+$", "", gsub('"', "", x))
  to_num <- function(x) {
    x <- trim(x)
    x[x %in% c("", "NA", "NaN", "null", "NULL", "#N/A")] <- NA
    suppressWarnings(as.numeric(x))
  }
  guess_delim <- function(path) {
    ln <- readLines(path, n = 5L, warn = FALSE)
    ln <- ln[nzchar(ln)][1L]
    if (lengths(regmatches(ln, gregexpr("\t", ln, fixed = TRUE))) >=
        lengths(regmatches(ln, gregexpr(",", ln, fixed = TRUE)))) "\t" else ","
  }
  read_char_matrix <- function(path, delim = guess_delim(path)) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    sp <- strsplit(lines, delim, fixed = TRUE)
    nc <- max(lengths(sp))
    m <- do.call(rbind, lapply(sp, function(x) {
      length(x) <- nc
      x[is.na(x)] <- ""
      x
    }))
    m
  }
}

norm_name <- function(x) {
  x <- tolower(trim(x))
  gsub("[^a-z0-9]+", "_", x)
}

## Normalized names of non-sample metadata columns
.META_NAMES <- c(
  # MCR feature matrices
  "numbering_for_matrix", "numbering_for_unique_peaks", "uniquepeak",
  "representative_peak_selected", "marked_peaks", "iis_i_a_or_u", "i_a_m_or_p",
  "final_annotation", "ri", "win", "iamp", "tag", "lineno", "query", "known",
  "cas", "preferred", "action", "peak_nr", "peak", "scan_nr", "ret_umin",
  "ri_", "m_z", "mass", "mass_scan", "deconv", "annotation",
  "compound", "check", "mstag", "anplus", "parent", "child", "rt_min", "rt_s",
  "retention_time", "retention_time_min", "inchi", "smiles",
  "formula", "database", "name", "metabolites", "synonym", "biological_name",
  "num", "no", "index", "exp_level_norm", "expert_concentration",
  "expert_curated", "comment", "list", "box", "sample_type",
  "concentration", "conc", "design_name", "design", "operator_name",
  "batch_name", "unified_class", "level", "peak_intensity", "fold", "ave",
  "significant", "retention", "matched_ion", "mrm", "abs_x",
  "log_lower", "log_upper", "kinds", "rt_allin", "x", "cx", "pvalue",
  "adj_pvalue", "sample_id", "sample", "peak_name", "name_", "markerlynx_xs")

.HEADER_ANCHORS <- c(
  "win", "ri", "final annotation", "numbering for matrix", "iis_i_a_or_u",
  "ret(umin)", "peak nr", "annotation", "metabolites", "synonym", "name",
  "m_z", "m/z", "mass", "retention", "retention time")

#' Locate the header row in a character matrix (usually row 1; banner rows
#' for MarkerLynx exports). Returns 1-based row index.
find_header_row <- function(m, max_scan = 12L) {
  nr <- min(nrow(m), max_scan)
  for (i in seq_len(nr)) {
    row <- tolower(trim(m[i, ]))
    if (sum(nzchar(row)) < 3L) next
    if (any(row %in% .HEADER_ANCHORS) ||
        any(grepl("^ret\\(|^numbering|^peak nr|^m/z", row))) return(i)
  }
  for (i in seq_len(nr)) if (sum(nzchar(trim(m[i, ]))) >= 3L) return(i)
  1L
}

read_matrix_file <- function(path) {
  if (grepl("\\.xlsx?$|\\.xls$", path, ignore.case = TRUE)) {
    ex <- read_excel_matrix(path)
    list(path = path, hdr = ex$header, body = ex$body)
  } else {
    m <- read_char_matrix(path)
    if (nrow(m) < 2L) stop("matrix has no data rows: ", path)
    hr <- find_header_row(m)
    hdr <- trim(m[hr, ])
    body <- m[seq(hr + 1L, nrow(m)), , drop = FALSE]
    list(path = path, hdr = hdr, body = body)
  }
}

#' Classify matrix columns; returns list with indices/vectors used downstream.
classify_columns <- function(hdr, body) {
  nn <- norm_name(hdr)
  probe <- body[seq_len(min(nrow(body), 60L)), , drop = FALSE]

  num_frac <- vapply(seq_len(ncol(body)), function(j) {
    v <- to_num(probe[, j])
    frac <- mean(!is.na(v))
    ifelse(is.nan(frac), 0, frac)
  }, numeric(1))

  is_meta <- nn %in% .META_NAMES |
    grepl("^fdata|^numbering|^marked|^iamp|^synon|^opt_memory", nn)

  first_meta <- which(is_meta)
  if (!length(first_meta)) stop("no metadata columns recognized")
  last_meta <- max(first_meta)

  pick1 <- function(cands) {
    for (cd in cands) {
      w <- which(nn == cd)
      if (length(w)) return(w[1L])
    }
    NA_integer_
  }

  id_col    <- pick1(c("numbering_for_matrix", "numbering_for_unique_peaks",
                       "uniquepeak", "peak_nr", "peak", "metabolites", "name",
                       "synonym"))
  grade_col <- pick1(c("iis_i_a_or_u"))
  iamp_col  <- pick1(c("iamp", "i_a_m_or_p"))
  name_col  <- pick1(c("final_annotation", "synonym", "annotation", "compound",
                       "name", "metabolites", "biological_name", "peak_name"))
  win_col   <- pick1(c("win"))
  ri_col    <- pick1(c("ri", "ri_"))
  rt_col    <- pick1(c("ret_umin", "rt_min", "rt_s", "retention_time",
                       "retention_time_min", "rt"))
  mz_col    <- pick1(c("m_z", "mz", "mass"))

  ## RT unit from the header token
  rt_unit <- NA_character_
  if (!is.na(rt_col)) {
    t <- nn[rt_col]
    rt_unit <- if (grepl("umin|u|min", t)) "min" else
      if (grepl("rt_s|_s|second", t)) "s" else "min"
  }

  sample_idx <- which(num_frac >= 0.5 & seq_along(nn) > last_meta)
  sample_idx <- sample_idx[nzchar(nn[sample_idx])]
  if (!length(sample_idx)) stop("no sample columns found")

  family <- if (!is.na(win_col)) "gcms_mcr" else
            if (!is.na(name_col)) "named" else "unknown"

  list(id_col = id_col, grade_col = grade_col, iamp_col = iamp_col,
       name_col = name_col, win_col = win_col, ri_col = ri_col, rt_col = rt_col,
       rt_unit = rt_unit, mz_col = mz_col, sample_idx = sample_idx,
       nn = nn, num_frac = num_frac, family = family, hdr = hdr,
       is_meta = is_meta, last_meta = last_meta)
}

#' Read matrix -> normalized list ready for mzTab-M building.
#' Retention times are returned in seconds.
read_matrix <- function(path) {
  x <- read_matrix_file(path)
  cls <- classify_columns(x$hdr, x$body)
  b <- x$body
  n <- nrow(b)

  get_col <- function(i) if (is.na(i)) rep("", n) else trim(b[, i])
  feature_id <- get_col(cls$id_col)
  feature_id <- ifelse(nzchar(feature_id), feature_id, paste0("feature_", seq_len(n)))
  feature_id <- make.unique(feature_id)

  grade <- get_col(cls$grade_col)
  iamp  <- get_col(cls$iamp_col)
  name  <- gsub('^"|"$', "", get_col(cls$name_col))
  win_raw <- get_col(cls$win_col)
  win <- ifelse(grepl("Win[0-9]+_", win_raw),
                sub("^.*?(Win[0-9]+_[A-Za-z0-9]+).*$", "\\1", win_raw),
                win_raw)
  ri    <- to_num(get_col(cls$ri_col))
  mz    <- to_num(get_col(cls$mz_col))
  rt_raw <- to_num(get_col(cls$rt_col))

  rt_s <- rep(NA_real_, n)
  if (!is.na(cls$rt_col)) {
    rt_s <- switch(cls$rt_unit,
                   min = rt_raw * 60,
                   s   = rt_raw,
                   rt_raw * 60)
  } else if (!is.na(cls$win_col)) {
    ## MCR matrices: RT in seconds is embedded as the leading number of "Win"
    hit <- grepl("^\\s*[0-9]", win_raw)
    rt_s_raw <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", win_raw)))
    rt_s <- ifelse(hit & !is.na(rt_s_raw), rt_s_raw, NA_real_)
  }

  abund <- apply(b[, cls$sample_idx, drop = FALSE], 2, to_num)
  storage.mode(abund) <- "double"
  colnames(abund) <- make.unique(trim(x$hdr[cls$sample_idx]))
  rownames(abund) <- feature_id

  list(
    path = path,
    family = cls$family,
    feature_id = feature_id,
    feature_name = name,
    grade = grade,
    iamp = iamp,
    win_token = win,
    ri = ri,
    rt_s = rt_s,
    mz = mz,
    abund = abund,
    sample_names = colnames(abund),
    cls = cls
  )
}

#' Read phenodata and align it to the matrix sample columns.
read_phenodata <- function(path, sample_names) {
  if (grepl("\\.xlsx?$|\\.xls$", path, ignore.case = TRUE)) {
    ex <- read_excel_matrix(path)
    hdr <- ex$header
    body <- ex$body
  } else {
    m <- read_char_matrix(path)
    hr <- find_header_row(m)
    hdr <- trim(m[hr, ])
    body <- m[seq(hr + 1L, nrow(m)), , drop = FALSE]
  }
  if (!ncol(body) || !length(hdr)) stop("empty phenodata: ", path)
  body <- body[, seq_along(hdr), drop = FALSE]
  colnames(body) <- make.unique(ifelse(nzchar(trim(hdr)),
                                       trim(hdr), paste0("V", seq_along(hdr))))
  df <- as.data.frame(body, stringsAsFactors = FALSE)

  res <- vapply(seq_len(ncol(df)), function(j) {
    v <- trim(df[[j]])
    bonus <- if (norm_name(colnames(df)[j]) %in%
        c("gcid", "sample_id", "sample", "samplename", "sample_name", "id",
          "filename", "sample id", "ms_sp", "msp")) 3 else 0
    sum(v %in% sample_names) + bonus
  }, numeric(1))
  idc <- which.max(res)
  id_col <- colnames(df)[idc]
  ids <- trim(df[[id_col]])
  hit_rate <- mean(sample_names %in% ids)
  if (hit_rate < 0.3) stop("phenodata sample-ID column mismatch: ", path)

  jr <- vapply(seq_len(ncol(df)), function(j) {
    v <- trim(df[[j]])
    v <- v[nzchar(v)]
    u <- unique(v)
    length(u) > 0 && length(u) <= 3 && all(nchar(u) <= 40) &&
      !all(v %in% c("null", "NA")) && j != idc
  }, logical(1))
  groups <- colnames(df)[jr]

  list(df = df, id_col = id_col, ids = ids,
       matched = mean(ids %in% sample_names), groups = groups)
}
