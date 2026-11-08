input_dir <- Sys.getenv("HONEYBEE_INPUT_DIR", "inputs")
output_dir <- Sys.getenv("HONEYBEE_FIGURE_DIR", "figures")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
if (length(args) > 1L || !mode %in% c("all", "counts", "ideogram", "ideogram-display")) {
  stop("Usage: Rscript 03_repertoire_genome.R [all|counts|ideogram|ideogram-display]", call. = FALSE)
}

if (mode %in% c("all", "counts")) {
  path <- file.path(input_dir, "repertoire", "receptor_counts.tsv")
  if (!file.exists(path)) {
    stop("Missing final repertoire counts: ", path, ". Run with 'ideogram-display' for Figure S1C.",
         call. = FALSE)
  }
  counts <- read.delim(path, colClasses = "character", check.names = FALSE)
  species <- c("A.mellifera", "D.melanogaster")
  families <- c("IR", "OR")
  needed <- c("species", "family", "total_annotated", "olfactory_detected")
  if (!all(needed %in% names(counts)) || anyDuplicated(names(counts)) ||
      nrow(counts) != 4L || !setequal(paste(counts$species, counts$family),
                                     as.vector(outer(species, families, paste)))) {
    stop("Expected one count row per A.mellifera/D.melanogaster and IR/OR.", call. = FALSE)
  }
  for (nm in needed[3:4]) counts[[nm]] <- suppressWarnings(as.numeric(counts[[nm]]))
  values <- as.matrix(counts[needed[3:4]])
  if (any(!is.finite(values) | values < 0 | values != floor(values)) ||
      any(counts$olfactory_detected > counts$total_annotated)) {
    stop("Repertoire counts must be nonnegative integers and detected cannot exceed total.",
         call. = FALSE)
  }
  draw_counts <- function() {
    op <- par(mfrow = c(2, 1), mar = c(4.6, 3.4, 1.0, 4.9), mgp = c(1.65, 0.38, 0),
              xpd = NA, tcl = -0.22, cex = 0.80)
    on.exit(par(op))
    for (i in seq_along(families)) {
      d <- counts[counts$family == families[i], ]
      d <- d[match(species, d$species), ]
      ink <- if (i == 1) c("#FFD6AB", "#E85712") else c("#CFDBF2", "#5982C2")
      ymax <- max(c(70, 160)[i], d$total_annotated)
      plot.new()
      plot.window(c(0.48, 2.55), c(0, ymax), xaxs = "i", yaxs = "i")
      rect(1:2 - 0.23, 0, 1:2 + 0.23, d$total_annotated, col = ink[1], border = NA)
      rect(1:2 - 0.18, 0, 1:2 + 0.18, d$olfactory_detected, col = ink[2], border = NA)
      axis(1, at = c(0.5, 2.55), labels = FALSE, lwd.ticks = 0)
      axis(2, at = seq(0, ymax, by = c(20, 40)[i]), las = 1)
      text(1:2, -0.06 * ymax, species, srt = 28, adj = 1, font = 3)
      mtext("Count", side = 2, line = 2.25)
      legend("topright", inset = c(-0.87, 0), bty = "n", cex = 0.83, x.intersp = 0.55,
             fill = ink, border = NA, legend = c("Total annotated", "Olfactory detected"))
    }
  }
}

if (mode != "counts") {
  bands <- read.delim("data/repertoire/ideogram_bands.tsv", colClasses = "character")
  loci <- read.delim("data/repertoire/ir_loci.tsv", stringsAsFactors = FALSE)
  hex <- strsplit(bands$colors, ",", fixed = TRUE)
  numeric_columns <- c("chr", "end_bp", "x_pt", "top_pt", "bottom_pt")
  bands[numeric_columns] <- lapply(bands[numeric_columns], as.numeric)
  if (!identical(bands$chr, as.numeric(1:16)) ||
      !all(lengths(hex) == ceiling(bands$end_bp / 200000)) ||
      !all(grepl("^[0-9A-F]{6}$", unlist(hex))) ||
      any(!is.finite(as.matrix(bands[numeric_columns]))) ||
      any(bands$top_pt >= bands$bottom_pt) || nrow(loci) != 11L ||
      anyDuplicated(loci$feature) || any(!loci$chr %in% bands$chr) ||
      any(loci$start_bp < 1 | loci$start_bp > loci$end_bp |
          loci$end_bp > bands$end_bp[match(loci$chr, bands$chr)])) {
    stop("Invalid chromosome band or study IR locus table.", call. = FALSE)
  }
  colors <- grDevices::colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100)
  numeric_density <- mode != "ideogram-display"
  if (numeric_density) {
    path <- file.path(input_dir, "repertoire", "transcript_density.tsv")
    if (!file.exists(path)) {
      stop("Missing 200-kb transcript-start counts: ", path,
           ". Run with 'ideogram-display' to use the verified display bands.", call. = FALSE)
    }
    density <- read.delim(path, colClasses = "numeric", check.names = FALSE)
    needed <- c("chr", "start", "end", "count")
    if (!all(needed %in% names(density)) || anyDuplicated(names(density)) ||
        any(!is.finite(as.matrix(density[needed]))) ||
        any(as.matrix(density[needed]) != floor(as.matrix(density[needed]))) ||
        any(density$count < 0) || !setequal(density$chr, bands$chr) ||
        max(density$count) < 1) {
      stop("Invalid density table; expected integer chr, start, end, count.", call. = FALSE)
    }
    hex <- lapply(seq_len(nrow(bands)), function(i) {
      d <- density[density$chr == i, ]
      d <- d[order(d$start), ]
      starts <- seq(0, bands$end_bp[i], by = 200000)
      if (!identical(d$start, starts) ||
          !identical(d$end, pmin(starts + 200000, bands$end_bp[i]))) {
        stop("Missing or unexpected 200-kb windows on chromosome ", i, call. = FALSE)
      }
      sub("^#", "", colors[1 + round(99 * d$count / max(density$count))])
    })
    bands$top_pt <- bands$bottom_pt - 160 * bands$end_bp / max(bands$end_bp)
    i <- match(loci$chr, bands$chr)
    loci$y_pt <- bands$top_pt[i] + (bands$bottom_pt[i] - bands$top_pt[i]) *
      ((loci$start_bp + loci$end_bp) / 2) / bands$end_bp[i]
  }
  draw_ideogram <- function() {
    op <- par(mar = rep(0, 4))
    on.exit(par(op))
    plot.new()
    plot.window(c(12, 366), c(352, 175), xaxs = "i", yaxs = "i", asp = 1)
    width <- 5.62
    radius <- 1.65
    for (i in seq_len(nrow(bands))) {
      d <- bands[i, ]
      edges <- if (numeric_density) {
        d$top_pt + c(seq(0, d$end_bp, by = 200000), d$end_bp) /
          d$end_bp * (d$bottom_pt - d$top_pt)
      } else {
        seq(d$top_pt, d$bottom_pt, length.out = length(hex[[i]]) + 1L)
      }
      inset <- function(y) {
        cap <- pmin(y - d$top_pt, d$bottom_pt - y)
        radius - sqrt(pmax(0, radius^2 - pmax(0, radius - cap)^2))
      }
      for (j in seq_along(hex[[i]])) {
        yy <- seq(edges[j], edges[j + 1L], length.out = 7)
        dx <- inset(yy)
        polygon(c(d$x_pt + dx, rev(d$x_pt + width - dx)), c(yy, rev(yy)),
                col = paste0("#", hex[[i]][j]), border = NA)
      }
      yy <- sort(c(seq(d$top_pt, d$top_pt + radius, length.out = 12),
                   seq(d$bottom_pt - radius, d$bottom_pt, length.out = 12)))
      dx <- inset(yy)
      polygon(c(d$x_pt + dx, rev(d$x_pt + width - dx)), c(yy, rev(yy)),
              border = "#7A8588", lwd = 0.8)
      text(d$x_pt + width / 2, d$bottom_pt + 7, i, srt = 90, cex = 0.66)
    }
    rect(loci$x_pt, loci$y_pt - 1.30, loci$x_pt + 2.8, loci$y_pt + 1.30,
         col = "#DEA99D", border = NA)
    x <- seq(290, 342, length.out = length(colors) + 1L)
    rect(head(x, -1), 200.5, tail(x, -1), 207.5, col = colors, border = NA)
    text(c(290, 342), 212.2, c("Low", "High"), cex = 0.7)
    rect(297, 221, 303.5, 227.5, col = "#DEA99D", border = NA)
    text(308, 224.3, "Ir genes", pos = 4, offset = 0, cex = 0.73)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(plot, file, width, height) {
  grDevices::pdf(file.path(output_dir, file), width, height, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  plot()
}
if (mode %in% c("all", "counts")) {
  # Figure S1A
  save_plot(draw_counts, "FigureS1A_receptor_counts.pdf", 3.6, 3.7)
}
if (mode != "counts") {
  # Figure S1C
  save_plot(draw_ideogram, "FigureS1C.pdf", 6.4, 3.2)
}
