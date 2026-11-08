input_dir <- Sys.getenv("HONEYBEE_INPUT_DIR", "inputs")
output_dir <- Sys.getenv("HONEYBEE_FIGURE_DIR", "figures")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
if (length(args) > 1L || !mode %in% c("all", "fig2", "s4")) {
  stop("Usage: Rscript 02_locus_tracks.R [all|fig2|s4]", call. = FALSE)
}

read_input <- function(file, columns) {
  path <- file.path(input_dir, "osn", file)
  if (!file.exists(path)) stop("Missing ", path, ". See INPUTS.md.", call. = FALSE)
  x <- read.delim(path, colClasses = "character", check.names = FALSE)
  if (!all(columns %in% names(x)) || anyDuplicated(names(x))) {
    stop("Expected columns in ", path, ": ", paste(columns, collapse = ", "), call. = FALSE)
  }
  x
}

as_numbers <- function(x, columns) {
  for (column in columns) {
    x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
    if (any(!is.finite(x[[column]]))) stop("Invalid numeric column: ", column, call. = FALSE)
  }
  x
}

main_genes <- c("Orco", "Ir76b", "Ir75g")
other_genes <- c("Ir25a", "Ir8a", "Ir76a", "Ir75e", "Ir75f")
specific <- list(Ir76a = c("0", "5"), Ir75e = "0", Ir75f = "5")
used <- switch(mode, fig2 = main_genes, s4 = other_genes, c(main_genes, other_genes))
tracks <- read_input("locus_coverage.tsv", c("gene", "group", "chr", "position_bp", "signal",
                                           "y_min", "y_max"))
tracks <- as_numbers(tracks[tracks$gene %in% used, ],
                     c("position_bp", "signal", "y_min", "y_max"))
exons <- read_input("locus_exons.tsv", c("gene", "start_bp", "end_bp", "strand"))
exons <- as_numbers(exons[exons$gene %in% used, ], c("start_bp", "end_bp"))
if (!setequal(tracks$gene, used) || !setequal(exons$gene, used) ||
    any(tracks$position_bp < 1 | tracks$signal < 0 | tracks$y_min < 0 |
        tracks$y_min >= tracks$y_max) ||
    anyDuplicated(tracks[c("gene", "group", "position_bp")]) ||
    any(exons$start_bp < 1 | exons$start_bp > exons$end_bp) ||
    any(!exons$strand %in% c("+", "-"))) {
  stop("Check the coverage values, displayed limits and exon coordinates.", call. = FALSE)
}
if (mode != "fig2") {
  cells <- read_input("osn_cells.tsv", c("cell", "cluster", names(specific)))
  cells <- as_numbers(cells, names(specific))
  if (nrow(cells) != 4175L || anyDuplicated(cells$cell) ||
      !setequal(cells$cluster, as.character(0:59)) ||
      sum(cells$cluster == "0") != 55L || sum(cells$cluster == "5") != 90L ||
      any(as.matrix(cells[names(specific)]) < 0)) {
    stop("Expected the final RNA table with 4,175 nuclei and clusters 0-59.", call. = FALSE)
  }
}
palette <- read.delim("data/osn/cluster_palette.tsv", colClasses = "character", comment.char = "")
cluster_ink <- setNames(palette$color, palette$cluster)
group_ink <- c(positive = "#dd6e42", negative = "#e8dab2")
cluster_order <- c(sort(as.character(1:59), method = "radix"), "0")

library(ggplot2)

coverage_plot <- function(gene, label_groups = TRUE) {
  d <- tracks[tracks$gene == gene, ]
  e <- exons[exons$gene == gene, ]
  is_main <- gene %in% main_genes
  groups <- if (is_main) cluster_order else names(group_ink)
  limits <- unique(d[c("group", "y_min", "y_max")])
  x_limits <- range(d$position_bp)
  if (!setequal(d$group, groups) || nrow(limits) != length(groups) ||
      min(table(d$group)) < 2 || length(unique(d$chr)) != 1 ||
      length(unique(e$strand)) != 1 || diff(x_limits) <= 0 ||
      any(e$end_bp < x_limits[1] | e$start_bp > x_limits[2])) {
    stop("Unexpected groups, axis limits or annotation for ", gene, call. = FALSE)
  }
  d <- d[order(d$group, d$position_bp), ]
  d$baseline <- length(groups) + 1 - match(d$group, groups)
  d$height <- d$baseline + .86 * (pmin(pmax(d$signal, d$y_min), d$y_max) - d$y_min) /
    (d$y_max - d$y_min)
  limits$y <- length(groups) + 1 - match(limits$group, groups)
  limits$label <- paste0("[", limits$y_min, " - ", limits$y_max, "]")
  group_labels <- if (is_main) groups else paste0(gene, c("+ OSNs", "- OSNs"))
  start <- max(x_limits[1], min(e$start_bp))
  end <- min(x_limits[2], max(e$end_bp))
  arrow <- grid::arrow(length = grid::unit(1.5, "mm"), type = "closed",
                       ends = if (e$strand[1] == "+") "last" else "first")
  p <- ggplot(d, aes(position_bp)) +
    geom_ribbon(aes(ymin = baseline, ymax = height, fill = group, group = group), colour = NA) +
    scale_fill_manual(values = if (is_main) cluster_ink else group_ink, guide = "none") +
    annotate("segment", x = start, xend = end, y = .40, yend = .40, arrow = arrow,
             colour = "#536476", size = .24) +
    geom_rect(data = e, aes(xmin = start_bp, xmax = end_bp, ymin = .24, ymax = .56),
              inherit.aes = FALSE, fill = "#536476") +
    scale_x_continuous(labels = if (is_main) function(x) format(x/1000, trim = TRUE)
                        else function(x) format(x, scientific = FALSE, trim = TRUE),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = rev(seq_along(groups)) + .24, labels = group_labels,
                       expand = c(0, 0)) +
    coord_cartesian(xlim = x_limits, ylim = c(0, length(groups) + 1), clip = "on") +
    labs(title = gene, y = if (label_groups) "Normalized signal" else NULL,
         x = paste0(d$chr[1], if (is_main) " position (kb)" else " position (bp)")) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "italic", hjust = .5),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = if (label_groups) element_text(size = if (is_main) 5 else 7)
                        else element_blank())
  if (!is_main) {
    p <- p + geom_text(data = limits, aes(x = x_limits[2], y = y + .69, label = label),
                       inherit.aes = FALSE, hjust = 1.05, size = 2.1)
  } else {
    p <- p + labs(subtitle = paste(unique(limits$label), collapse = ", "))
  }
  p
}

expression_plot <- function(gene) {
  positive <- cells$cluster %in% specific[[gene]]
  d <- data.frame(value = cells[[gene]], group = factor(
    ifelse(positive, "positive", "negative"), levels = c("negative", "positive")))
  set.seed(1)
  take <- unlist(lapply(split(seq_len(nrow(d)), d$group), function(i) {
    if (length(i) > 800L) sample(i, 800L) else i
  }), use.names = FALSE)
  labels <- c(paste0(gene, "- OSNs (n = ", sum(!positive), ")"),
              paste0(gene, "+ OSNs (n = ", sum(positive), ")"))
  ggplot(d, aes(value, group, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .85, size = .35) +
    geom_boxplot(width = .16, fill = "white", outlier.shape = NA, size = .25) +
    geom_point(data = d[take, ], position = position_jitter(width = 0, height = .05, seed = 1),
               size = .18, alpha = .12) +
    scale_fill_manual(values = group_ink, guide = "none") +
    scale_y_discrete(labels = labels) +
    labs(x = paste("Exp. level of", gene), y = NULL) +
    theme_classic(base_size = 8)
}

save_panels <- function(panels, file, width, height, vertical = FALSE) {
  grDevices::pdf(file.path(output_dir, file), width, height, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  n <- length(panels)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = if (vertical) grid::grid.layout(n, 1)
                                   else grid::grid.layout(1, n)))
  for (i in seq_len(n)) {
    vp <- grid::viewport(layout.pos.row = if (vertical) i else 1,
                         layout.pos.col = if (vertical) 1 else i)
    print(panels[[i]], newpage = FALSE, vp = vp)
  }
  grid::popViewport()
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (mode != "s4") {
  panels <- lapply(seq_along(main_genes), function(i) coverage_plot(main_genes[i], i == 1))
  # Figure 2A
  save_panels(panels, "Figure2A_accessibility_tracks.pdf", 5.3, 10.4)
}
if (mode != "fig2") {
  # Figure S4B
  save_panels(list(coverage_plot("Ir25a")), "FigureS4B_Ir25a_accessibility.pdf", 3.3, 2.6)
  # Figure S4C
  save_panels(list(coverage_plot("Ir8a")), "FigureS4C_Ir8a_accessibility.pdf", 3.3, 2.6)
  # Figure S4D
  save_panels(list(coverage_plot("Ir76a"), expression_plot("Ir76a")),
              "FigureS4D_Ir76a_accessibility_expression.pdf", 3.3, 4, TRUE)
  # Figure S4E
  save_panels(list(coverage_plot("Ir75e"), expression_plot("Ir75e")),
              "FigureS4E_Ir75e_accessibility_expression.pdf", 3.3, 4, TRUE)
  # Figure S4F
  save_panels(list(coverage_plot("Ir75f"), expression_plot("Ir75f")),
              "FigureS4F_Ir75f_accessibility_expression.pdf", 3.3, 4, TRUE)
}
