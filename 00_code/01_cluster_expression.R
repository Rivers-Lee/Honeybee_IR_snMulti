input_dir <- Sys.getenv("HONEYBEE_INPUT_DIR", "inputs")
output_dir <- Sys.getenv("HONEYBEE_FIGURE_DIR", "figures")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
if (length(args) > 1L || !mode %in% c("all", "evidence", "violins")) {
  stop("Usage: Rscript 01_cluster_expression.R [all|evidence|violins]", call. = FALSE)
}

genes <- c("Orco", "Ir76b", "Ir75g", "Ir25a", "Ir8a", "Ir75f", "Ir76a", "Ir75e")
calls <- read.delim("data/osn/figure4b_calls.tsv", colClasses = "character")
if (nrow(calls) != 60L || !setequal(calls$cluster, as.character(0:59)) ||
    anyNA(calls$calls) || any(!grepl("^[LSA.]{8}$", calls$calls))) {
  stop("Expected 60 unique rows of eight final Figure 4B calls.", call. = FALSE)
}
codes <- do.call(rbind, strsplit(calls$calls, "", fixed = TRUE))
dimnames(codes) <- list(calls$cluster, genes)
coreceptors <- c("Orco", "Ir76b", "Ir25a", "Ir8a")
signature <- apply(codes[, coreceptors] != ".", 1, function(active) {
  paste(coreceptors[active], collapse = "+")
})
class <- ifelse(signature == "Ir76b+Ir25a", "C1",
         ifelse(signature %in% c("Orco+Ir76b+Ir25a", "Orco+Ir76b+Ir8a"), "C2",
         ifelse(signature == "Orco+Ir76b+Ir25a+Ir8a", "C3", NA_character_)))
classes <- c("C1", "C2", "C3")
counts <- table(factor(class, levels = classes))
if (anyNA(class) || !identical(as.integer(counts), c(1L, 19L, 40L))) {
  stop("Final integrated co-receptor calls do not yield the published C1/C2/C3 counts.",
       call. = FALSE)
}

if (mode != "evidence") {
  path <- file.path(input_dir, "osn", "osn_cells.tsv")
  if (!file.exists(path)) {
    stop("Missing final RNA/counts table: ", path, ". Run with 'evidence' for Figure 4B/C.",
         call. = FALSE)
  }
  cells <- read.delim(path, check.names = FALSE, colClasses = "character")
  targets <- genes[1:5]
  values <- c(targets, paste0(targets, "_UMI"))
  if (!all(c("cell", "cluster", values) %in% names(cells)) || anyDuplicated(names(cells)) ||
      nrow(cells) != 4175L || anyNA(cells$cell) || any(!nzchar(cells$cell)) ||
      anyDuplicated(cells$cell) || !setequal(cells$cluster, calls$cluster) ||
      sum(cells$cluster == "0") != 55L || sum(cells$cluster == "5") != 90L) {
    stop("Expected final 4,175 nuclei with cluster IDs 0-59 and RNA/UMI columns.", call. = FALSE)
  }
  for (nm in values) {
    x <- suppressWarnings(as.numeric(cells[[nm]]))
    if (any(!is.finite(x) | x < 0) || (grepl("_UMI$", nm) && any(x != floor(x)))) {
      stop("Invalid RNA data or counts in ", nm, call. = FALSE)
    }
    cells[[nm]] <- x
  }
}

library(ggplot2)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (mode != "violins") {
  palette <- c(C1 = "#224C4B", C2 = "#42846E", C3 = "#85A68B")
  d <- expand.grid(cluster = calls$cluster, gene = genes, stringsAsFactors = FALSE)
  d$call <- as.vector(codes)
  d$x <- match(d$gene, genes)
  d$y <- 61 - match(d$cluster, calls$cluster)
  d$ink <- ifelse(d$gene == "Orco", "#2A5C8B", "#C74B39")
  strip <- data.frame(y = 60:1, colour = unname(palette[class]))
  p <- ggplot(d, aes(x, y)) +
    geom_tile(fill = "white", colour = "grey83", size = 0.15) +
    geom_vline(xintercept = c(0.5, 3.5, 5.5, 8.5), size = 0.28) +
    geom_hline(yintercept = c(0.5, 60.5, 61.8), size = 0.28) +
    geom_point(data = d[d$call == "L", ], aes(colour = ink), size = 1.9) +
    geom_point(data = d[d$call == "S", ], aes(colour = ink), size = 0.95) +
    geom_point(data = d[d$call == "A", ], aes(colour = ink), shape = 21,
               fill = "white", size = 0.85, stroke = 0.18) +
    geom_point(data = d[d$call == ".", ], shape = 4, colour = "grey50", size = 0.72,
               stroke = 0.20) +
    geom_tile(data = strip, aes(x = 8.70, y = y, fill = colour),
              width = 0.28, height = 1, inherit.aes = FALSE) +
    scale_colour_identity() + scale_fill_identity() +
    annotate("text", x = c(2, 4.5, 7), y = 61.15,
             label = c("Universal", "Broad", "Specific"), size = 2.45) +
    scale_x_continuous(breaks = 1:8, labels = genes, limits = c(0.47, 8.9),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 60:1, labels = calls$cluster, limits = c(0.4, 61.85),
                       expand = c(0, 0)) +
    labs(caption = paste("Large: mean > 0.1 and detection > 10%; small: mean > 0 and detection < 10%",
                         "Open: accessible only; cross: below criteria", sep = "\n")) +
    theme_classic(base_size = 8) +
    theme(axis.title = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1,
                                                                  face = "italic"),
          axis.text.y = element_text(size = 6), plot.caption = element_text(hjust = 0, size = 6))

  # Figure 4B
  ggsave(file.path(output_dir, "Figure4B.pdf"), p, width = 3.9, height = 8.7,
         device = "pdf", useDingbats = FALSE)

  assignments <- data.frame(cluster = calls$cluster, class = class, complement = signature)
  summary <- aggregate(list(n_clusters = rep(1L, nrow(assignments))),
                       assignments[c("class", "complement")], sum)
  # Figure 4C
  write.table(assignments, file.path(output_dir, "Figure4C_cluster_assignments.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Figure 4C
  write.table(summary, file.path(output_dir, "Figure4C_complements.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

if (mode != "evidence") {
  colors <- read.delim("data/osn/cluster_palette.tsv", colClasses = "character", comment.char = "")
  cluster_palette <- setNames(colors$color, colors$cluster)
  common_order <- names(cluster_palette)
  broad_order <- list(
    Ir25a = c(47,24,37,22,16,45,21,7,26,42,3,17,25,23,28,15,38,31,14,36,27,32,10,4,
              18,29,6,40,46,8,44,39,52,2,30,48,20,41,35,51,56,43,53,59,54,19,11,12,
              58,13,55,0,34,9,1,33,49,5,50,57),
    Ir8a = c(31,21,25,4,15,41,12,3,10,50,22,33,43,59,6,46,58,5,40,44,24,57,37,42,
             48,23,20,17,36,13,26,49,18,16,19,27,34,35,39,53,9,0,1,11,14,2,38,28,
             45,29,47,30,51,32,52,54,55,56,7,8))
  ticks <- list(Orco = c(1, 3, 10, 30), Ir76b = c(1, 3, 5), Ir75g = 1:3,
                Ir25a = c(1, 3, 5), Ir8a = 1:3)
  star_y <- c(Ir76b = 1.55, Ir75g = 1.55, Ir25a = 1.85, Ir8a = 1.6)
  violin <- function(gene, order, labels = TRUE) {
    umi_sum <- tapply(cells[[paste0(gene, "_UMI")]], cells$cluster, sum)
    cutoff <- if (gene == "Orco") 4 else 1
    raw <- cells[[gene]]
    raw[umi_sum[cells$cluster] < cutoff] <- 0
    d <- data.frame(x = factor(cells$cluster, levels = order), value = raw + 1)
    star_clusters <- rownames(codes)[codes[, gene] == "A"]
    p <- ggplot(d, aes(x, value, fill = x)) +
      geom_violin(scale = "width", trim = TRUE, size = 0.17) +
      scale_fill_manual(values = cluster_palette, guide = "none") +
      scale_y_log10(breaks = ticks[[gene]]) +
      labs(x = NULL, y = NULL, title = gene) +
      theme_classic(base_size = 8) +
      theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 8),
            axis.text.x = if (labels) element_text(angle = 45, hjust = 1, size = 6) else element_blank(),
            axis.ticks.x = if (labels) element_line() else element_blank(),
            plot.margin = margin(1, 4, 1, 3))
    if (length(star_clusters)) {
      p <- p + geom_point(data = data.frame(x = factor(star_clusters, levels = order),
                                           y = unname(star_y[gene])),
                         aes(x, y), inherit.aes = FALSE, shape = 8, size = 1.2,
                         stroke = 0.32, colour = "#E56C23")
    }
    p
  }
  save_stack <- function(plots, file, height) {
    grDevices::pdf(file.path(output_dir, file), width = 10, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
    grid::grid.newpage()
    grid::grid.text("Expression Level", x = 0.015, rot = 90, gp = grid::gpar(fontsize = 8))
    grid::pushViewport(grid::viewport(x = 0.52, width = 0.96,
                                     layout = grid::grid.layout(length(plots), 1)))
    for (i in seq_along(plots)) {
      print(plots[[i]], newpage = FALSE, vp = grid::viewport(layout.pos.row = i))
    }
    grid::popViewport()
  }
  # Figure 1F
  save_stack(list(violin("Orco", common_order, FALSE), violin("Ir76b", common_order, FALSE),
                  violin("Ir75g", common_order)), "Figure1F.pdf", 3.5)
  # Figure S4A
  save_stack(lapply(names(broad_order), function(g) violin(g, as.character(broad_order[[g]]))),
             "FigureS4A.pdf", 2.7)
}
