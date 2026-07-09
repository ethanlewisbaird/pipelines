#!/usr/bin/env Rscript
# dump_seurat.R
# Stream every matrix and metadata field out of a Seurat .rds object to disk,
# without densifying sparse matrices or building a giant in-memory copy.
#
# Usage: Rscript dump_seurat.R <input.rds> <out_dir>
#
# Output: <out_dir>/manifest.json  + raw binary component files per matrix.

suppressWarnings(suppressMessages({
  library(Matrix)
  library(jsonlite)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript dump_seurat.R <input.rds> <out_dir> [assay_name]")
}
rds_path   <- args[[1]]
out_dir    <- args[[2]]
force_assay <- if (length(args) >= 3) args[[3]] else NULL
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf(...), "\n", file = stderr())

msg("Reading %s ...", rds_path)
obj <- readRDS(rds_path)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

manifest <- list(
  object_class = class(obj)[1],
  layers = list(), obsm = list(), varm = list(),
  obsp = list(), obs = list(), var = list(), uns = list()
)

# Write a sparse matrix as raw CSC components (no densify, no copy beyond cast).
write_sparse <- function(m, tag) {
  m <- as(m, "CsparseMatrix")  # dgCMatrix (CSC)
  writeBin(as.double(m@x),    file.path(out_dir, paste0(tag, ".x.f64")), size = 8, endian = "little")
  writeBin(as.integer(m@i),   file.path(out_dir, paste0(tag, ".i.i32")), size = 4, endian = "little")
  writeBin(as.integer(m@p),   file.path(out_dir, paste0(tag, ".p.i32")), size = 4, endian = "little")
  list(format = "csc", tag = tag, shape = dim(m), nnz = length(m@x))
}

# Write a dense matrix in column-major (R native) order.
write_dense <- function(m, tag) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  writeBin(as.double(m), file.path(out_dir, paste0(tag, ".d.f64")), size = 8, endian = "little")
  list(format = "dense", tag = tag, shape = dim(m), order = "F")
}

write_matrix <- function(m, tag) {
  if (is.null(m) || prod(dim(m)) == 0) return(NULL)
  if (inherits(m, "sparseMatrix") || inherits(m, "Matrix")) {
    out <- tryCatch(write_sparse(m, tag), error = function(e) {
      msg("  sparse write failed for %s (%s); writing dense", tag, conditionMessage(e))
      write_dense(m, tag)
    })
  } else {
    out <- write_dense(m, tag)
  }
  out
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

has_slot <- function(x, name) isTRUE(tryCatch(.hasSlot(x, name), error = function(e) FALSE))

# ---------------------------------------------------------------------------
# Resolve default assay
# ---------------------------------------------------------------------------

assay_names <- tryCatch(names(obj@assays), error = function(e) character(0))
if (length(assay_names) == 0) stop("No assays found on object.")

default_assay <- tryCatch({
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    SeuratObject::DefaultAssay(obj)
  } else if (has_slot(obj, "active.assay")) {
    obj@active.assay
  } else assay_names[1]
}, error = function(e) assay_names[1])
if (is.null(default_assay) || !default_assay %in% assay_names) default_assay <- assay_names[1]

if (!is.null(force_assay)) {
  if (!force_assay %in% assay_names) stop(sprintf("Assay '%s' not found. Available: %s", force_assay, paste(assay_names, collapse=", ")))
  msg("Overriding default assay '%s' -> '%s'", default_assay, force_assay)
  default_assay <- force_assay
}

assay <- obj@assays[[default_assay]]
manifest$main_assay <- default_assay
msg("Default assay: %s  (classes: %s)", default_assay, paste(class(assay), collapse = ","))

# ---------------------------------------------------------------------------
# Extract expression slots — handles Seurat v4 (@counts/@data/@scale.data)
# and v5 (@layers list). For v5 split layers, join via name prefix.
# ---------------------------------------------------------------------------

get_v5_layers <- function(assay) {
  if (!has_slot(assay, "layers")) return(NULL)
  L <- assay@layers
  if (is.null(L) || length(L) == 0) return(NULL)
  L
}

collect_slots <- list()  # canonical name -> matrix

v5 <- get_v5_layers(assay)
if (!is.null(v5)) {
  msg("Seurat v5 layers detected: %s", paste(names(v5), collapse = ", "))
  # Group split layers (e.g. counts.1, counts.2) by their base name.
  base_of <- function(nm) sub("\\..*$", "", nm)
  bases <- unique(vapply(names(v5), base_of, character(1)))
  cells_target <- tryCatch(colnames(obj), error = function(e) NULL)
  for (b in bases) {
    parts <- names(v5)[vapply(names(v5), function(nm) base_of(nm) == b, logical(1))]
    if (length(parts) == 1) {
      collect_slots[[b]] <- v5[[parts]]
    } else {
      msg("  joining %d split layers for '%s'", length(parts), b)
      mats <- lapply(parts, function(p) as(v5[[p]], "CsparseMatrix"))
      joined <- tryCatch(do.call(cbind, mats), error = function(e) {
        msg("  cbind failed for %s (%s); using first part only", b, conditionMessage(e))
        mats[[1]]
      })
      # Reorder columns to global cell order when possible.
      if (!is.null(cells_target) && !is.null(colnames(joined)) &&
          all(cells_target %in% colnames(joined))) {
        joined <- joined[, cells_target, drop = FALSE]
      }
      collect_slots[[b]] <- joined
    }
  }
} else {
  for (sn in c("counts", "data", "scale.data")) {
    m <- tryCatch(slot(assay, sn), error = function(e) NULL)
    if (!is.null(m) && length(m) > 0 && prod(dim(m)) > 0) collect_slots[[sn]] <- m
  }
}

if (length(collect_slots) == 0) stop("No expression matrix found in default assay.")

# Pick a reference slot for dimnames (prefer data, then counts, then first).
ref_name <- if ("data" %in% names(collect_slots)) "data" else
            if ("counts" %in% names(collect_slots)) "counts" else
            names(collect_slots)[1]
ref <- collect_slots[[ref_name]]
manifest$var_names <- rownames(ref) %||% NULL
manifest$obs_names <- colnames(ref) %||% NULL
n_genes <- nrow(ref); n_cells <- ncol(ref)
msg("Reference slot '%s': %d genes x %d cells", ref_name, n_genes, n_cells)

for (sn in names(collect_slots)) {
  tag <- paste0("layer__", gsub("[^A-Za-z0-9]", "_", sn))
  meta <- write_matrix(collect_slots[[sn]], tag)
  if (!is.null(meta)) {
    manifest$layers[[gsub("[^A-Za-z0-9]", "_", sn)]] <- meta
    msg("  wrote layer %s [%s %dx%d]", sn, meta$format, meta$shape[1], meta$shape[2])
  }
  collect_slots[[sn]] <- NULL  # free
  gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Reductions: embeddings (cells x dims) -> obsm; loadings (genes x dims) -> varm
# ---------------------------------------------------------------------------

red_names <- tryCatch(names(obj@reductions), error = function(e) character(0))
for (rname in red_names) {
  red <- obj@reductions[[rname]]
  emb <- tryCatch(red@cell.embeddings, error = function(e) NULL)
  if (!is.null(emb) && length(emb) > 0) {
    key <- paste0("X_", tolower(rname))
    meta <- write_matrix(emb, paste0("obsm__", gsub("[^A-Za-z0-9]", "_", key)))
    if (!is.null(meta)) { manifest$obsm[[key]] <- meta; msg("  wrote obsm %s", key) }
  }
  load <- tryCatch(red@feature.loadings, error = function(e) NULL)
  if (!is.null(load) && length(load) > 0 && prod(dim(load)) > 0) {
    key <- paste0(toupper(rname), "s")
    meta <- write_matrix(load, paste0("varm__", gsub("[^A-Za-z0-9]", "_", key)))
    if (!is.null(meta)) { manifest$varm[[key]] <- meta; msg("  wrote varm %s", key) }
  }
  gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Graphs / neighbor matrices (cells x cells) -> obsp
# ---------------------------------------------------------------------------

graph_names <- tryCatch(names(obj@graphs), error = function(e) character(0))
for (gname in graph_names) {
  g <- obj@graphs[[gname]]
  if (!is.null(g) && length(g) > 0 && all(dim(g) == n_cells)) {
    key <- gsub("[^A-Za-z0-9]", "_", gname)
    meta <- write_matrix(g, paste0("obsp__", key))
    if (!is.null(meta)) { manifest$obsp[[key]] <- meta; msg("  wrote obsp %s", key) }
  }
  gc(verbose = FALSE)
}

# Also capture @neighbors if present (Seurat stores some graphs here).
nn_names <- tryCatch(names(obj@neighbors), error = function(e) character(0))
for (nn in nn_names) {
  g <- tryCatch(obj@neighbors[[nn]], error = function(e) NULL)
  if (!is.null(g) && inherits(g, "Matrix") && all(dim(g) == n_cells)) {
    key <- gsub("[^A-Za-z0-9]", "_", paste0("nn_", nn))
    meta <- write_matrix(g, paste0("obsp__", key))
    if (!is.null(meta)) { manifest$obsp[[key]] <- meta; msg("  wrote obsp %s", key) }
  }
}

# ---------------------------------------------------------------------------
# obs (cell metadata)
# ---------------------------------------------------------------------------

md <- tryCatch(obj@meta.data, error = function(e) NULL)
if (!is.null(md) && ncol(md) > 0) {
  for (k in colnames(md)) {
    col <- md[[k]]
    if (is.factor(col)) {
      manifest$obs[[k]] <- list(
        type = "categorical",
        values = as.integer(col),     # 1-based codes, NA preserved
        levels = levels(col)
      )
    } else if (is.logical(col)) {
      manifest$obs[[k]] <- list(type = "bool", values = col)
    } else {
      manifest$obs[[k]] <- list(type = "array", values = col)
    }
  }
  msg("Captured %d obs columns", ncol(md))
}

# ---------------------------------------------------------------------------
# var (feature metadata)
# ---------------------------------------------------------------------------

mf <- tryCatch(assay@meta.features, error = function(e) NULL)
if (!is.null(mf) && ncol(mf) > 0) {
  for (k in colnames(mf)) {
    col <- mf[[k]]
    if (is.factor(col)) {
      manifest$var[[k]] <- list(type = "categorical",
                                values = as.integer(col), levels = levels(col))
    } else if (is.logical(col)) {
      manifest$var[[k]] <- list(type = "bool", values = col)
    } else {
      manifest$var[[k]] <- list(type = "array", values = col)
    }
  }
  msg("Captured %d var columns", ncol(mf))
}

# Variable features list, if any.
vf <- tryCatch(VariableFeatures(obj), error = function(e) NULL)
if (is.null(vf)) vf <- tryCatch(assay@var.features, error = function(e) NULL)
if (!is.null(vf) && length(vf) > 0) manifest$uns$variable_features <- as.character(vf)

# ---------------------------------------------------------------------------
# uns: project name, command history
# ---------------------------------------------------------------------------

manifest$uns$project_name <- tryCatch(as.character(obj@project.name), error = function(e) NULL)
cmds <- tryCatch(names(obj@commands), error = function(e) NULL)
if (!is.null(cmds) && length(cmds) > 0) manifest$uns$seurat_commands <- as.character(cmds)

# ---------------------------------------------------------------------------
# Write manifest
# ---------------------------------------------------------------------------

write_json(manifest, file.path(out_dir, "manifest.json"),
           auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
msg("Manifest written. Done.")
