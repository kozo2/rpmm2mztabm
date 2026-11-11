## io.R - low-level file reading: delimiter detection, matrix scanning, excel helpers
##
## All RPMM data files are text (.txt/.tsv/.csv) or Excel (.xlsx/.xls).
## Delimiters cannot be trusted from extensions (e.g. RPMM0016/gcms_INTEGRATED.csv
## is tab-delimited), so detection is content-based.

#' @return "\t" or "," based on the first non-empty line.
guess_delim <- function(path) {
  ln <- readLines(path, n = 5L, warn = FALSE)
  ln <- ln[nzchar(ln)][1L]
  if (is.na(ln)) stop("empty file: ", path)
  if (lengths(regmatches(ln, gregexpr("\t", ln, fixed = TRUE))) >=
      lengths(regmatches(ln, gregexpr(",", ln, fixed = TRUE)))) "\t" else ","
}

#' Read a delimited file into a character matrix (no name repair, no factors),
#' with consistent row lengths (ragged rows padded).
read_char_matrix <- function(path, delim = guess_delim(path), header_rows = 1L) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\x1a?$", lines)]           # drop empty / DOS-EOF lines
  if (!length(lines)) stop("no content: ", path)
  sp <- strsplit(lines, delim, fixed = TRUE)
  nc <- max(lengths(sp))
  m <- do.call(rbind, lapply(sp, function(x) {
    length(x) <- nc
    x[is.na(x)] <- ""
    x
  }))
  m
}

#' Trim whitespace; also strip the padding seen in some files (e.g. Kusano
#' "R_01_1           ").
trim <- function(x) gsub("^\\s+|\\s+$", "", gsub('"', "", x))

#' Convert character -> numeric, "" / "NA" / "null" -> NA.
#' Not global: it only trims whitespace and surrounding quotes.
to_num <- function(x) {
  x <- trim(x)
  x[x %in% c("", "NA", "NaN", "null", "NULL", "#N/A")] <- NA
  suppressWarnings(as.numeric(x))
}

#' Excel reading with header row discovery.
#' Many sheets carry banner rows and meta rows before the real header
#' (e.g. RPMM0013: 4 meta rows then the matrix header). We scan for the row
#' containing a known anchor token ('Win' or 'RI' or 'final annotation' or
#' 'Numbering for matrix') and use it as the header.
read_excel_matrix <- function(path, sheet = 1L) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("readxl is required for ", path)
  sh <- readxl::excel_sheets(path)
  sn <- if (is.character(sheet)) sheet else sh[sheet]
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = sn, col_names = FALSE, col_types = "text"))
  m <- as.matrix(raw)
  storage.mode(m) <- "character"
  m[is.na(m)] <- ""
  anchors <- c("Win", "RI", "final annotation", "Numbering for matrix",
               "IIS_I_A_or_U", "Metabolite", "retention", "m/z")
  hdr_row <- NA_integer_
  for (i in seq_len(min(nrow(m), 25L))) {
    row <- trim(m[i, ])
    if (sum(nzchar(row)) < 3L) next
    if (any(tolower(row) %in% tolower(anchors)) ||
        any(grepl("^[Rr]et\\b|Ret\\(|Numbering", row))) {
      hdr_row <- i
      break
    }
  }
  if (is.na(hdr_row)) {
    ## fall back: first row with >= 3 non-empty cells
    for (i in seq_len(min(nrow(m), 25L))) if (sum(nzchar(trim(m[i, ]))) >= 3L) {
      hdr_row <- i
      break
    }
  }
  if (is.na(hdr_row)) stop("cannot locate header in ", path, " sheet ", sn)
  list(header = trim(m[hdr_row, ]),
       body   = m[seq(hdr_row + 1L, nrow(m)), , drop = FALSE],
       pre_header = if (hdr_row > 1L) m[seq_len(hdr_row - 1L), , drop = FALSE] else NULL,
       sheet = sn, sheets = sh)
}

#' Column-wise numeric-frame check: how many cells look numeric.
numeric_fraction <- function(col) {
  v <- to_num(col)
  mean(!is.na(v))
}

#' Exact monoisotopic masses for the simple formula parser (named-metabolite LUT).
.proton <- 1.0072764666370876    # m/z of [M+H]+ uses monoisotopic H+ (1.00727647)
.element_mass <- c(
  C = 12.0, H = 1.00782503223, N = 14.00307400443, O = 15.99491461957,
  S = 31.9720711744, P = 30.97376199842, Si = 27.97692653503,
  Cl = 34.968852682, Br = 78.9183376, K = 38.9637064864, Na = 22.9897692820,
  Fe = 55.93493633, Ca = 39.962590863, Mg = 23.985041697, Zn = 63.92914201,
  Mn = 54.93804391, Cu = 62.92959772, Se = 79.9165218, F = 18.99840320677)

#' Parse a molecular formula like "C40H56O2" and return the [M+H]+ mass.
mass_plus_H <- function(formula) {
  if (!nzchar(formula) || is.na(formula)) return(NA_real_)
  toks <- regmatches(formula, gregexpr("[A-Z][a-z]?[0-9]*", formula))[[1L]]
  if (paste(toks, collapse = "") != formula) return(NA_real_)
  el <- sub("[0-9]*$", "", toks)
  n <- suppressWarnings(as.numeric(sub("^[A-Z][a-z]?", "", toks)))
  n[is.na(n)] <- 1
  ok <- el %in% names(.element_mass)
  if (!all(ok)) return(NA_real_)
  sum(.element_mass[el] * n) + .proton
}
