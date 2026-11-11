## parse_msp.R - read NIST MSP annotation files (*.txt starting with "Name:")
##
## RPMM MSP files hold one block per deconvoluted GC-MS feature. The block
## key is the feature "Win###_C###" token which also appears in the matrix
## "Win" column, so matrices and MSP blocks join on that token.
##
## Example block:
##   Name: 230.63  Win001_C01
##   Synon:Annotated:
##   Synon:Classified:
##   Synon:RI:997
##   Synon:RT(s):230.6
##   Synon:MCR-Reference: Win001_C01
##   Synon:Protocol:Freeze dry sample
##   Synon:Project:130224_Chiba_U_tomato2011
##   Synon:Species:tomato Leaf, Fruit and stem
##   DB#: 1
##   Num Peaks:6
##    73  999; 79  864; 103   46; 145  979;
##    162   64; 173   58;
##
## Alkane calibration files ("ALKANE_SERIE.txt") use the same container with
## "Name: C10" etc. - they carry RI anchors for provenance only.

#' Parse one MSP-like file -> data.frame (one row per block).
parse_msp <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^Name:", lines)
  if (!length(starts)) return(NULL)
  ends <- c(starts[-1L] - 1L, length(lines))
  blocks <- lapply(seq_along(starts), function(i) lines[starts[i]:ends[i]])

  syn <- function(blk, key) {
    hit <- grep(paste0("^Synon:(", key, "):"), blk, value = TRUE, ignore.case = TRUE)
    if (!length(hit)) return(NA_character_)
    v <- trim(sub("^[^:]*:[^:]*:", "", hit))   # drop "Synon:<key>:"
    ifelse(nzchar(v) | grepl("annotated|classified", key), v, NA)
  }

  one <- function(blk) {
    name <- sub("^Name:\\s*", "", blk[1L])
    name <- trim(name)
    npi <- grep("^Num Peaks:", blk)
    npeaks <- if (length(npi)) as.numeric(sub("^Num Peaks:\\s*", "", blk[npi[1L]])) else NA_real_
    pk <- if (length(npi)) blk[(npi[1L] + 1L):length(blk)] else character()
    pk <- pk[grepl("\\d\\s+[0-9]", pk)]
    toks <- trimws(unlist(strsplit(paste(pk, collapse = " "), ";", fixed = TRUE)))
    toks <- toks[nzchar(toks)]
    bp_mz <- bp_int <- NA_real_
    if (length(toks)) {
      mm <- do.call(rbind, strsplit(toks, "\\s+"))
      if (ncol(mm) >= 2L) {
        mz <- suppressWarnings(as.numeric(mm[, 1L]))
        int <- suppressWarnings(as.numeric(mm[, 2L]))
        k <- which.max(int)
        if (length(k)) {
          bp_mz <- mz[k]
          bp_int <- int[k]
        }
      }
    }
    ri  <- suppressWarnings(as.numeric(syn(blk, "RI")))
    rts <- suppressWarnings(as.numeric(syn(blk, "RT\\(s\\)|RT")))
    data.frame(
      name = name,
      win_token = sub("^.*?(Win[0-9]+_[A-Za-z0-9]+).*$", "\\1", name),
      ri = if (length(ri)) ri[1L] else NA_real_,
      rt_s = if (length(rts)) rts[1L] else NA_real_,
      base_peak_mz = bp_mz,
      base_peak_int = bp_int,
      n_peaks = npeaks,
      db = trim(sub("^DB#:\\s*", "", blk[grep("^DB#:", blk)[1L] %||% NA])),
      annotated = syn(blk, "Annotated"),
      classified = syn(blk, "Classified"),
      pi = syn(blk, "PI"),
      species = syn(blk, "Species"),
      protocol = syn(blk, "Protocol"),
      project = syn(blk, "Project"),
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, lapply(blocks, one))
  df$win_token[!grepl("^Win[0-9]+_", df$win_token)] <- NA
  rownames(df) <- NULL
  df
}

`%||%` <- function(a, b) if (is.null(a) || !length(a) || is.na(a)) b else a

#' Index MSP blocks by Win token (first occurrence wins when duplicates exist).
msp_win_index <- function(msp_df) {
  if (is.null(msp_df) || !nrow(msp_df)) return(list())
  df <- msp_df[!is.na(msp_df$win_token), , drop = FALSE]
  df <- df[!duplicated(df$win_token), , drop = FALSE]
  split(df, df$win_token)
}

#' Pick m/z for feature rows given matrix Win/RI and an MSP index.
#' Priority: exact Win token match -> nearest MSP block by RI -> NA.
solve_feature_mz <- function(win_token, ri, msp_df) {
  n <- length(win_token)
  mz <- rep(NA_real_, n)
  src <- rep(NA_character_, n)
  idx <- msp_win_index(msp_df)

  for (i in seq_len(n)) {
    tok <- win_token[i]
    if (!is.na(tok) && !is.null(idx[[tok]])) {
      row <- idx[[tok]]
      if (!is.na(row$base_peak_mz)) {
        mz[i] <- row$base_peak_mz
        src[i] <- "base peak m/z of deconvoluted EI spectrum (Win token match)"
        next
      }
    }
  }
  miss <- which(is.na(mz))
  if (length(miss) && !is.null(msp_df) && nrow(msp_df) &&
      any(!is.na(msp_df$base_peak_mz))) {
    good <- which(!is.na(msp_df$ri) & !is.na(msp_df$base_peak_mz))
    for (i in miss) {
      if (is.na(ri[i]) || !length(good)) next
      k <- good[which.min(abs(msp_df$ri[good] - ri[i]))]
      mz[i] <- msp_df$base_peak_mz[k]
      src[i] <- "base peak m/z of deconvoluted EI spectrum (nearest-RI match)"
    }
  }
  list(mz = mz, source = src)
}
