## inventory.R - classify every file in investigation/downloads/data and pair
## matrices to phenodata + NIST MSP annotation files.
##
## Fails-and-logs: anything the classifiers cannot resolve is recorded in a
## curation table rather than guessed. Per-project species CV terms are only
## assigned from a curated map of strings actually observed in MSP
## "Synon:Species:" lines and phenodata columns; unknown species leave colData
## $species unset and land in the curation list.

source(file.path("investigation/mztabm/lib", "io.R"))
source(file.path("investigation/mztabm/lib", "read_matrix.R"))
source(file.path("investigation/mztabm/lib", "parse_msp.R"))

## ---- curated species map --------------------------------------------------
## Only taxonomy IDs I am confident in; every raw string is logged with the
## mapping decision so reviewers can correct it in projects.tsv.
.species_map <- c(
  "arabidopsis"           = "[NCBITaxon, NCBITaxon:3702, Arabidopsis thaliana, ]",
  "arabi"                 = "[NCBITaxon, NCBITaxon:3702, Arabidopsis thaliana, ]",
  "tomato"                = "[NCBITaxon, NCBITaxon:4081, Solanum lycopersicum, ]",
  "rice"                  = "[NCBITaxon, NCBITaxon:4530, Oryza sativa, ]",
  "ine"                   = "[NCBITaxon, NCBITaxon:4530, Oryza sativa, ]",
  "lettuce"               = "[NCBITaxon, NCBITaxon:4236, Lactuca sativa, ]",
  "sunnylettuce"          = "[NCBITaxon, NCBITaxon:4236, Lactuca sativa, ]",
  "cassava"               = "[NCBITaxon, NCBITaxon:3983, Manihot esculenta, ]",
  "chickpea"              = "[NCBITaxon, NCBITaxon:3827, Cicer arietinum, ]",
  "soybean"               = "[NCBITaxon, NCBITaxon:3847, Glycine max, ]",
  "fermented soybean"     = "[NCBITaxon, NCBITaxon:3847, Glycine max, ]",
  "soy"                   = "[NCBITaxon, NCBITaxon:3847, Glycine max, ]",
  "natane"                = "[NCBITaxon, NCBITaxon:3708, Brassica napus, ]",
  "kabu"                  = "[NCBITaxon, NCBITaxon:3711, Brassica rapa, ]",
  "brassica rapa"         = "[NCBITaxon, NCBITaxon:3711, Brassica rapa, ]",
  "brassica napus"        = "[NCBITaxon, NCBITaxon:3708, Brassica napus, ]")

## Strings seen in the data whose taxonomy ID is not certain enough to emit:
## they are logged to the curation list and species is left unset.
.species_unresolved <- c(
  "aspen", "kanabikisou", "allium plants", "allium species", "fistulosum",
  "chinense", "tuberosum", "spring onion", "green spring onion", "scallion",
  "leek", "rakkyo", "red spring onion", "oriental garlic", "white spring onion",
  "ine_r", "ine_root", "ine leaf", "stds", "arabi", "seeds", "fruit and stem",
  "tomato leaf, fruit and stem", "arabi and rice", "cassava and rice",
  "soy nodule leaf and root", "chickpea (leaf,root and nodule)",
  "fermented soybean", "arabi , seeds", "lettuce leaf")

#' Resolve a raw species string to a CV parameter or NA (logged, not guessed).
resolve_species <- function(raw) {
  v <- tolower(trim(raw))
  v <- gsub("\\s+", " ", v)
  v <- sub("\\s*\\(.*$", "", v)                       # drop parentheticals
  v <- sub("\\s+leaf.*$|\\s+fruit.*$|\\s+root.*$|\\s+seed.*$|\\s+stem.*$|\\s+nodule.*$",
           "", v)
  v <- trim(v)
  if (!nzchar(v)) return(NA_character_)
  if (v %in% names(.species_map)) return(unname(.species_map[v]))
  ## multi-species strings ("tomato leaf, fruit and stem" after reduction)
  for (k in names(.species_map)) if (startsWith(v, k)) return(unname(.species_map[k]))
  NA_character_
}

#' Guess scan polarity; GC-MS (EI) measures positive ions throughout.
guess_polarity <- function(relpath) {
  p <- tolower(relpath)
  if (grepl("neg|negative", p)) return("negative")
  if (grepl("pos|positive", p)) return("positive")
  "positive"
}

.is_limma <- function(fn) grepl("limma|diff_res|_res\\.", fn, ignore.case = TRUE)
.is_rdata <- function(fn) grepl("\\.rdata$|\\.rda$", fn, ignore.case = TRUE)
.is_alkane <- function(fn) grepl("alkane", fn, ignore.case = TRUE)

#' Classify one file: matrix | phenodata | msp | alkane | companion | skip
classify_file <- function(path) {
  fn <- basename(path)
  if (fn == ".DS_Store") return("skip")
  if (.is_rdata(fn)) return("skip")
  if (.is_limma(fn)) return("companion")           # differential results: out of scope
  if (grepl("^phenodata|^pheno", fn, ignore.case = TRUE)) return("phenodata")
  ext <- tolower(tools::file_ext(fn))
  if (!ext %in% c("txt", "csv", "tsv", "xlsx", "xls")) return("companion")

  if (!grepl("\\.xlsx?$|\\.xls$", fn, ignore.case = TRUE)) {
    first <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                      error = function(e) "")
    if (length(first) && grepl("^Name:", first)) {
      return(if (.is_alkane(fn)) "alkane" else "msp")
    }
  }
  ## everything left is a candidate quantification matrix; inventory calls
  ## read_header_info() to confirm and familises it
  "matrix_candidate"
}

#' Cheap header probe: family + sample-ness (no full parse). Returns NA family
#' when the file cannot be a quant matrix.
probe_matrix <- function(path) {
  tryCatch({
    x <- read_matrix_file(path)
    if (ncol(x$body) < 2L || nrow(x$body) < 1L) return(NA_character_)
    cls <- classify_columns(x$hdr, x$body)
    cls$family
  }, error = function(e) NA_character_)
}

#' Pair each matrix with the phenodata + MSP to use for it.
#' Priority: same directory -> sub/sibling dirs of the same AnalysedDataset
#' -> anywhere in the project. Ambiguity (several candidates) is logged.
pick_companion <- function(matrix_path, kind, files_df) {
  md <- dirname(matrix_path)
  proj <- sub("^([^/]+)/.*", "\\1",
              sub("^investigation/downloads/data/", "", matrix_path))
  cand <- files_df$path[files_df$kind == kind & files_df$project == proj]
  if (!length(cand)) return(list(path = NA_character_, note = "none in project"))
  same_dir <- cand[dirname(cand) == md]
  if (length(same_dir) == 1L) return(list(path = same_dir, note = "same dir"))
  if (length(same_dir) > 1L) {
    ## prefer basename stem overlap with the matrix stem
    stem <- tolower(gsub("[^a-z0-9]+", "", tools::file_path_sans_ext(basename(matrix_path))))
    ov <- vapply(same_dir, function(p) {
      s <- tolower(gsub("[^a-z0-9]+", "", tools::file_path_sans_ext(basename(p))))
      nchar(s) > 0 && (grepl(substr(s, 1, min(8, nchar(s))), stem) ||
                         grepl(substr(stem, 1, min(8, nchar(stem))), s))
    }, logical(1))
    if (sum(ov) == 1L) return(list(path = same_dir[ov], note = "same dir, name match"))
    return(list(path = same_dir[1L],
                note = paste0("AMBIGUOUS same dir (", length(same_dir), "); took first")))
  }
  ## project-wide: closest directory
  if (length(cand) == 1L) return(list(path = cand, note = "project-unique"))
  depth <- nchar(gsub("[^/]", "", dirname(cand)))
  md_depth <- nchar(gsub("[^/]", "", md))
  j <- order(abs(depth - md_depth))
  list(path = cand[j[1L]],
       note = paste0("project-wide nearest of ", length(cand), " candidate(s)"))
}

#' Build the full conversion inventory.
#'
#' @param data_root  root of the RPMM data tree
#' @param projects   named character vector: RPMM ID -> title (from projects.json)
#' @return list(manifest, curation, files)
build_inventory <- function(data_root = "investigation/downloads/data",
                            projects = NULL) {
  all_paths <- list.files(data_root, recursive = TRUE, full.names = TRUE)
  all_paths <- all_paths[!dir.exists(all_paths)]
  rel <- sub(paste0("^", data_root, "/"), "", all_paths)

  kind <- vapply(seq_along(all_paths), function(i) {
    classify_file(all_paths[i])
  }, character(1))

  cat("inventory: classified", length(all_paths), "files\n")
  print(table(kind))

  ## confirm matrix candidates and discover their family
  family <- rep(NA_character_, length(all_paths))
  ex <- which(kind == "matrix_candidate")
  cat("probing", length(ex), "matrix candidates ...\n")
  for (i in ex) {
    family[i] <- probe_matrix(all_paths[i])
    if (is.na(family[i])) {
      kind[i] <- "companion"
    } else {
      kind[i] <- "matrix"
    }
  }

  files_df <- data.frame(
    path = all_paths,
    rel = rel,
    project = sub("/.*", "", rel),
    kind = kind,
    family = family,
    stringsAsFactors = FALSE)

  ## ---- species per project (from MSP Synon:Species, then phenodata) ------
  curation <- list()
  add_curation <- function(project, issue, detail) {
    curation[[length(curation) + 1L]] <<- data.frame(
      project = project, issue = issue, detail = detail, stringsAsFactors = FALSE)
  }

  species_cv <- setNames(rep(NA_character_, length(unique(files_df$project))),
                         sort(unique(files_df$project)))
  species_raw <- setNames(rep(NA_character_, length(species_cv)), names(species_cv))

  for (proj in names(species_cv)) {
    msps <- files_df$path[files_df$project == proj & files_df$kind == "msp"]
    raw <- character()
    for (mp in head(msps, 3L)) {
      sp <- tryCatch({
        lns <- readLines(mp, warn = FALSE)
        v <- sub("^Synon:Species:", "", grep("^Synon:Species:", lns, value = TRUE))
        unique(trim(v)[nzchar(trim(v))])
      }, error = function(e) character())
      raw <- c(raw, sp)
    }
    if (!length(raw)) {
      phs <- files_df$path[files_df$project == proj & files_df$kind == "phenodata"]
      for (ph in phs) {
        d <- tryCatch({
          m <- read_char_matrix(ph)
          hdr <- norm_name(m[1L, ])
          j <- which(grepl("species|binomial", hdr))
          if (!length(j)) character() else
            unique(trim(m[-1L, j, drop = FALSE][nzchar(m[-1L, j])]))
        }, error = function(e) character())
        raw <- c(raw, d)
      }
    }
    raw <- unique(raw[nzchar(raw)])
    if (!length(raw)) {
      add_curation(proj, "species_missing",
                   "no species string in MSP or phenodata; colData$species omitted")
      next
    }
    species_raw[proj] <- paste(raw, collapse = " | ")
    cvs <- unique(na.omit(vapply(raw, resolve_species, character(1))))
    if (length(cvs) == 1L) {
      species_cv[proj] <- cvs
    } else if (length(cvs) > 1L) {
      add_curation(proj, "species_conflict",
                   paste0("multiple species -> ", paste(raw, collapse = " | "),
                          "; resolved to: ", paste(cvs, collapse = " | ")))
      species_cv[proj] <- cvs[1L]     # per-file override if phenodata has per-sample species
    } else {
      add_curation(proj, "species_unresolved_taxon",
                   paste0("raw string '", paste(raw, collapse = " | "),
                          "' has no curated NCBITaxon mapping; colData$species omitted"))
    }
  }

  ## ---- project titles ---------------------------------------------------
  title_of <- function(proj) {
    if (is.null(projects)) return("")
    t <- projects[proj]
    if (length(t) && !is.na(t)) t else ""
  }

  ## ---- pair matrices ----------------------------------------------------
  mtx <- files_df[files_df$kind == "matrix", , drop = FALSE]
  rows <- vector("list", nrow(mtx))
  for (i in seq_len(nrow(mtx))) {
    mp <- mtx$path[i]
    proj <- mtx$project[i]
    ph <- pick_companion(mp, "phenodata", files_df)
    sp <- pick_companion(mp, "msp", files_df)

    notes <- c(
      if (grepl("^AMBIGUOUS", ph$note)) paste0("phenodata: ", ph$note),
      if (grepl("^AMBIGUOUS", sp$note)) paste0("msp: ", sp$note))
    for (nt in notes) add_curation(proj, "pairing_ambiguity",
                                   paste0(mtx$rel[i], ": ", nt))

    ## per-sample species from phenodata when a species column exists
    sp_cv <- species_cv[proj]

    rows[[i]] <- data.frame(
      project = proj,
      matrix_path = mp,
      matrix_rel = mtx$rel[i],
      family = mtx$family[i],
      phenodata_path = if (is.na(ph$path)) "" else ph$path,
      phenodata_note = ph$note,
      msp_path = if (is.na(sp$path)) "" else sp$path,
      msp_note = sp$note,
      polarity = guess_polarity(mtx$rel[i]),
      species_cv = ifelse(is.na(sp_cv), "", sp_cv),
      species_raw = species_raw[proj] %||% "",
      title = title_of(proj),
      stringsAsFactors = FALSE)
  }
  manifest <- do.call(rbind, rows)

  ## ---- named-metabolite LUT coverage probe ------------------------------
  for (i in which(manifest$family == "named")) {
    nm <- tryCatch(unique(read_matrix(manifest$matrix_path[i])$feature_name),
                   error = function(e) character())
    nm <- nm[nzchar(nm) & !is.na(nm)]
    add_curation(manifest$project[i], "named_metabolites",
                 paste0(manifest$matrix_rel[i], ": ", length(nm),
                        " distinct metabolite name(s): ",
                        paste(head(sort(nm), 40), collapse = "; ")))
  }

  list(manifest = manifest,
       curation = if (length(curation)) do.call(rbind, curation) else
         data.frame(project = character(), issue = character(),
                    detail = character(), stringsAsFactors = FALSE),
       files = files_df)
}
