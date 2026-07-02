#!/usr/bin/env Rscript
# Gene-level TWAS Manhattan plot via ggplot2 + ggrepel.
# Ported verbatim from EXP-GRN-TWAS/src/manhattan_ggplot.R (the renderer whose
# ggrepel labels do not overlap), with two additions:
#   * --cap   : ceiling applied to -log10(p) (default 50, as in the original;
#               keeps one extreme peak from compressing the rest).
#   * also writes a sibling .png next to the .pdf output.
# Called from utils/plot_manhattan.py via utils/rscript_wrapper.sh.
#
# Usage:
#   Rscript manhattan_ggplot.R \
#       --input  <csv>   # columns: gene_name, chrom, pos, p_twas, p_twas_adj
#       --output <pdf>
#       --title  <str>
#       --fdr    0.05
#       --color  "#E41A1C"   # colour for significant dots
#       --cap    50          # -log10(p) ceiling
#       --high-thresh 15     # label every gene with -log10(p) above this
#       --label-n 0          # plus the top-N significant genes (0 = none extra)
#       --width  18  --height 6

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(scales)
})

# -- Args ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) && i < length(args)) args[[i + 1]] else default
}

input_csv   <- get_arg("--input")
output_pdf  <- get_arg("--output",  "manhattan.pdf")
plot_title  <- get_arg("--title",   "Gene-level TWAS Manhattan")
fdr_thresh  <- as.numeric(get_arg("--fdr",   "0.05"))
sig_color   <- get_arg("--color",   "#E41A1C")
fig_width   <- as.numeric(get_arg("--width",  "18"))
fig_height  <- as.numeric(get_arg("--height", "6"))
label_n     <- as.integer(get_arg("--label-n",    "0"))
high_thresh <- as.numeric(get_arg("--high-thresh", "15"))
cap         <- as.numeric(get_arg("--cap", "50"))

if (is.null(input_csv)) stop("--input <csv> is required")

# -- Data ----------------------------------------------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE)

required_cols <- c("chrom", "pos", "p_twas", "p_twas_adj", "gene_name")
missing <- setdiff(required_cols, names(df))
if (length(missing)) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

# Keep autosomes only (1-22)
df <- df[!is.na(df$chrom) & !is.na(df$pos) & !is.na(df$p_twas), ]
df <- df[df$chrom >= 1 & df$chrom <= 22, ]
df <- df[order(df$chrom, df$pos), ]

# Build cumulative x-coordinates
chr_order <- 1:22
gap       <- 20e6
offsets   <- numeric(23)
cum <- 0.0
for (ch in chr_order) {
    offsets[ch] <- cum
    sub_pos <- df$pos[df$chrom == ch]
    if (length(sub_pos) == 0) next
    cum <- cum + as.numeric(max(sub_pos, na.rm = TRUE)) + as.numeric(gap)
}

df$x_plot <- df$pos + offsets[df$chrom]
df$neglog <- pmin(-log10(pmax(df$p_twas, 1e-300)), cap)
df$sig    <- df$p_twas_adj < fdr_thresh

# Chromosome midpoints (present chromosomes only)
chrs_present <- intersect(chr_order, unique(df$chrom))
chr_mids <- sapply(chrs_present, function(ch) {
    x <- df$x_plot[df$chrom == ch]; (min(x) + max(x)) / 2
})
names(chr_mids) <- as.character(chrs_present)
chr_labs <- as.character(chrs_present)

df$chr_col <- ifelse(df$chrom %% 2 == 1, "#4E79A7", "#A0CBE8")
thresh_y   <- -log10(fdr_thresh)

# -- Labels: every gene with -log10(p) > high_thresh, plus top label_n sig -----
sig_df   <- df[df$sig, ]
sig_df   <- sig_df[order(sig_df$p_twas), ]
high_df  <- df[df$neglog > high_thresh, ]
top_df   <- if (nrow(sig_df) > 0 && label_n > 0) head(sig_df, label_n) else data.frame()
label_df <- unique(rbind(high_df, top_df))

# -- Plot ----------------------------------------------------------------------
p <- ggplot(df, aes(x = x_plot, y = neglog)) +
    geom_point(data = df[!df$sig, ], aes(color = chr_col),
               size = 0.6, alpha = 0.5, shape = 16, show.legend = FALSE) +
    geom_point(data = df[df$sig, ], color = sig_color, size = 1.8, alpha = 0.95,
               shape = 21, fill = sig_color, stroke = 0.25, show.legend = FALSE) +
    geom_hline(yintercept = thresh_y, linetype = "dashed",
               color = "#CC0000", linewidth = 0.5, alpha = 0.7) +
    scale_color_identity() +
    scale_x_continuous(breaks = as.numeric(chr_mids), labels = chr_labs,
                       expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = plot_title,
         subtitle = sprintf("%d significant genes (FDR < %.2f)", sum(df$sig), fdr_thresh),
         x = "Chromosome", y = expression(-log[10](p))) +
    theme_classic(base_size = 12) +
    theme(
        plot.title    = element_text(size = 14, face = "bold",  hjust = 0),
        plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0),
        axis.text.x   = element_text(size = 11),
        axis.text.y   = element_text(size = 11),
        axis.title.x  = element_text(size = 12, margin = margin(t = 6)),
        axis.title.y  = element_text(size = 12, margin = margin(r = 6)),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin   = margin(8, 12, 8, 8)
    )

# ggrepel labels (adaptive size/padding for many labels)
nlab <- nrow(label_df)
lab_size   <- if (nlab > 600) 1.3 else if (nlab > 200) 1.9 else 3.5
lab_boxpad <- if (nlab > 600) 0.10 else 0.35
if (nlab > 0) {
    p <- p + geom_text_repel(
        data = label_df, aes(label = gene_name),
        size = lab_size, color = "#111111",
        segment.color = "#888888", segment.size = 0.2,
        box.padding = lab_boxpad, point.padding = 0.10,
        max.overlaps = Inf, force = 1.0,
        min.segment.length = 0.0, nudge_y = 0.4)
}

# -- Save pdf + png ------------------------------------------------------------
ggsave(output_pdf, plot = p, width = fig_width, height = fig_height,
       device = cairo_pdf, units = "in")
png_path <- sub("\\.pdf$", ".png", output_pdf)
ggsave(png_path, plot = p, width = fig_width, height = fig_height,
       dpi = 300, units = "in")
cat(sprintf("Saved %s and %s\n", output_pdf, png_path))
