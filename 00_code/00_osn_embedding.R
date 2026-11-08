input_dir <- Sys.getenv("HONEYBEE_INPUT_DIR", "inputs")
output_dir <- Sys.getenv("HONEYBEE_FIGURE_DIR", "figures")
path <- file.path(input_dir, "osn", "osn_cells.tsv")

genes <- c(Orco = "Or2", Ir76b = "LOC551704", Ir75g = "LOC102653640",
           Ir25a = "LOC412949", Ir8a = "LOC409777", Ir75f = "LOC102656088",
           Ir76a = "LOC725724", Ir75e = "LOC100576097")
args <- commandArgs(trailingOnly = TRUE)
export <- length(args) == 2L && identical(args[1], "export")
if (length(args) && !export) {
  stop("Usage: Rscript 00_osn_embedding.R [export path/to/final_OSN.rds]", call. = FALSE)
}
if (export) {
  if (!file.exists(args[2])) stop("Missing final OSN RDS: ", args[2], call. = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite ", path, call. = FALSE)
  obj <- readRDS(args[2])
  xy <- SeuratObject::Embeddings(obj, reduction = "tsne.rna")
  meta <- obj[[]]
  if (!"Annotation" %in% names(meta) || !setequal(rownames(xy), rownames(meta)) ||
      !all(c("rnaTSNE_1", "rnaTSNE_2") %in% colnames(xy))) {
    stop("Expected final Annotation metadata and tsne.rna coordinates.", call. = FALSE)
  }
  ids <- rownames(xy)
  cells <- data.frame(cell = ids, cluster = as.character(meta[ids, "Annotation"]),
                      rnaTSNE_1 = xy[, "rnaTSNE_1"], rnaTSNE_2 = xy[, "rnaTSNE_2"])
  cells$cluster[cells$cluster %in% "Orco-Neuron"] <- "0"
  for (layer in c("data", "counts")) {
    mat <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      SeuratObject::GetAssayData(obj, assay = "RNA", layer = layer)
    } else {
      SeuratObject::GetAssayData(obj, assay = "RNA", slot = layer)
    }
    if (!all(unname(genes) %in% rownames(mat)) || !setequal(ids, colnames(mat))) {
      stop("RNA/", layer, " lacks the expected matrix features or final cell IDs.", call. = FALSE)
    }
    for (g in names(genes)) {
      nm <- if (layer == "counts") paste0(g, "_UMI") else g
      cells[[nm]] <- as.numeric(mat[genes[[g]], ids])
      if (layer == "counts" && any(!is.finite(cells[[nm]]) |
                                  cells[[nm]] < 0 | cells[[nm]] != floor(cells[[nm]]))) {
        stop("Invalid original RNA counts for ", g, call. = FALSE)
      }
    }
  }
} else {
  if (!file.exists(path)) stop("Missing final t-SNE/RNA table: ", path, call. = FALSE)
  cells <- read.delim(path, check.names = FALSE, colClasses = "character")
}
needed <- c("cell", "cluster", "rnaTSNE_1", "rnaTSNE_2", names(genes))
if (!all(needed %in% names(cells)) || anyDuplicated(names(cells))) {
  stop("Expected unique columns: ", paste(needed, collapse = ", "), call. = FALSE)
}
if (nrow(cells) != 4175L || anyNA(cells$cell) || any(!nzchar(cells$cell)) ||
    anyDuplicated(cells$cell) || !setequal(cells$cluster, as.character(0:59)) ||
    sum(cells$cluster == "0") != 55L || sum(cells$cluster == "5") != 90L) {
  stop("Expected the final 4,175 nuclei, with unique cell IDs and clusters 0-59.", call. = FALSE)
}
for (nm in needed[-c(1:2)]) {
  cells[[nm]] <- suppressWarnings(as.numeric(cells[[nm]]))
  if (any(!is.finite(cells[[nm]]))) stop("Non-finite values in ", nm, call. = FALSE)
}
if (any(as.matrix(cells[names(genes)]) < 0) ||
    any(vapply(cells[needed[3:4]], function(x) diff(range(x)) == 0, logical(1)))) {
  stop("Expected nonnegative normalized RNA and two nonconstant t-SNE axes.", call. = FALSE)
}
if (export) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(cells, path, sep = "\t", row.names = FALSE, quote = TRUE)
  message("Exported existing t-SNE coordinates and RNA data/counts to ", path)
  quit(save = "no", status = 0)
}

library(ggplot2)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
xr <- range(cells$rnaTSNE_1)
yr <- range(cells$rnaTSNE_2)
dx <- diff(xr)
dy <- diff(yr)

feature_map <- function(gene, axes = FALSE, legend = TRUE, highlight = FALSE) {
  d <- cells[order(cells[[gene]]), ]
  d$display <- pmin(d[[gene]], 2)
  ink <- if (gene == "Orco") "#548CD4" else "#BF0000"
  p <- ggplot(d, aes(rnaTSNE_1, rnaTSNE_2, colour = display)) +
    geom_point(size = 0.42, stroke = 0) +
    scale_colour_gradient(low = "lightgrey", high = ink, limits = c(0, 2),
                          breaks = seq(0, 2, 0.5), name = "Expression") +
    coord_equal(xlim = xr + c(-0.075, 0.035) * dx, ylim = yr + c(-0.10, 0.035) * dy,
                expand = FALSE) +
    labs(title = gene) + theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "italic", size = 10, colour = ink),
          legend.position = if (legend) "right" else "none",
          legend.key.height = grid::unit(2.8, "mm"),
          legend.key.width = grid::unit(1.8, "mm"),
          plot.margin = margin(4, 4, 4, 4))
  if (axes) {
    x0 <- xr[1] - 0.035 * dx
    y0 <- yr[1] - 0.05 * dy
    a <- grid::arrow(length = grid::unit(1.8, "mm"), type = "closed")
    p <- p + annotate("segment", x = x0, xend = x0 + 0.28 * dx,
                      y = y0, yend = y0, arrow = a) +
      annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.30 * dy, arrow = a) +
      annotate("text", x = x0 + 0.14 * dx, y = y0 - 0.032 * dy, label = "rnaTSNE_1", size = 2.15) +
      annotate("text", x = x0 - 0.025 * dx, y = y0 + 0.15 * dy, label = "rnaTSNE_2",
               size = 2.15, angle = 90)
  }
  if (highlight) {
    cl0 <- cells[cells$cluster == "0", ]
    center <- c(mean(range(cl0$rnaTSNE_1)), mean(range(cl0$rnaTSNE_2)))
    radius <- max(sqrt((cl0$rnaTSNE_1 - center[1])^2 + (cl0$rnaTSNE_2 - center[2])^2)) +
      0.012 * max(dx, dy)
    theta <- seq(0, 2 * pi, length.out = 101)
    p <- p + annotate("path", x = center[1] + radius * cos(theta),
                      y = center[2] + radius * sin(theta), linetype = 2, size = 0.28) +
      annotate("text", x = center[1], y = center[2] + radius + 0.025 * dy,
               label = "0 (Orco- Neuron)", size = 2.2)
  }
  p
}

save_row <- function(plots, file, width, height) {
  grDevices::pdf(file.path(output_dir, file), width, height, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(plots))))
  for (i in seq_along(plots)) {
    print(plots[[i]], newpage = FALSE, vp = grid::viewport(layout.pos.col = i))
  }
  grid::popViewport()
}

# Figure 1C, Figure 1D, Figure 1E
save_row(list(feature_map("Orco", axes = TRUE, highlight = TRUE), feature_map("Ir76b"),
              feature_map("Ir75g")), "Figure1C-E.pdf", 8.5, 2.6)
# Figure 4A
save_row(lapply(1:5, function(i) feature_map(names(genes)[i + 3], i == 1, i == 5)),
         "Figure4A.pdf", 8.5, 1.9)
