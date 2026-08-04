###################################################################################################
##F1 Barcode Filtering Script (with barcode IDs, P1-like, P2-like, Hybrid)
###################################################################################################

# Load libraries
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)

# Read arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 4){
    stop("Usage: Rscript script.R [reads_count_file] [output_prefix] [p.call] [barcode_filter_file]")
}
counts <- as.character(args[1])            # input reads count file
name <- as.character(args[2])              # output prefix
p.call <- as.numeric(args[3])              # p.call (e.g., 0.9)
barcode_filter_file <- as.character(args[4]) # high_QC barcode file, obtained from Socrates analysis


###################################################################################################
## Step 1: Load read count data and apply barcode filter
###################################################################################################

# Read input count file
counts.m <- read.table(counts, header = TRUE, stringsAsFactors = FALSE)

# Read barcode filter file (use only first column)
filter.barcodes <- read.table(barcode_filter_file, header = TRUE, stringsAsFactors = FALSE)[[1]]

# Filter input count data by barcode intersection
counts.m <- counts.m[counts.m$barcode %in% filter.barcodes, ]

# Save filtered barcodes
barcodes <- counts.m$barcode

# Keep only reads columns
counts.m <- counts.m[, c("P1_reads", "P2_reads")]

# Ensure reads are numeric
counts.m$P1_reads <- as.numeric(counts.m$P1_reads)
counts.m$P2_reads <- as.numeric(counts.m$P2_reads)


###################################################################################################
## Step 2: Naive Bayes classifier
###################################################################################################

bayes <- function(data, p0 = 0.05, p1 = 0.05, pn = 0.9, p.call = 0.9) {
    data <- as.matrix(data)
    res <- apply(data, 1, function(x) {
        hap0 <- x[1]
        hap1 <- x[2]
        total <- hap0 + hap1
        p0_bc.loge <- dbinom(hap0, size = total, prob = p.call, log = TRUE)
        p1_bc.loge <- dbinom(hap1, size = total, prob = p.call, log = TRUE)
        p2_bc.loge <- dbinom(hap1, size = total, prob = 0.5, log = TRUE)

        p0_bc <- dbinom(hap0, size = total, prob = p.call)
        p1_bc <- dbinom(hap1, size = total, prob = p.call)
        p2_bc <- dbinom(hap1, size = total, prob = 0.5)

        p0_top <- p0_bc * p0
        p1_top <- p1_bc * p1
        p2_top <- p2_bc * pn
        total_prob <- p0_top + p1_top + p2_top
        prob_0 <- p0_top / total_prob
        prob_1 <- p1_top / total_prob
        prob_2 <- p2_top / total_prob

        return(list(P1.binom.loge = p0_bc.loge, P2.binom.loge = p1_bc.loge, Hybrid.binom.loge = p2_bc.loge,
                    P1.p = prob_0, P2.p = prob_1, Hybrid.p = prob_2))
    })
    res <- as.data.frame(do.call(rbind, res))
    return(res)
}

# Run Bayesian classification
res <- bayes(counts.m, p.call = p.call)
res <- as.data.frame(lapply(res, as.numeric))

###################################################################################################
## Step 3: Assign cells based on bayes test
###################################################################################################

# Assign based on max posterior
res$bayes.type <- apply(res[, c("P1.p", "P2.p", "Hybrid.p")], 1, function(x) { names(x)[which.max(x)] })
res$bayes.type <- gsub(".p", "", res$bayes.type)

# Attach counts and barcode
res$P1_reads <- counts.m$P1_reads
res$P2_reads <- counts.m$P2_reads
res$total <- res$P1_reads + res$P2_reads
res$max.p <- apply(res[, c("P1.p", "P2.p", "Hybrid.p")], 1, max)
res$barcode <- barcodes

# Move barcode to first column
res <- res[, c("barcode", setdiff(colnames(res), "barcode"))]

###################################################################################################
## Step 4: Filter perfect F1 cells
###################################################################################################

# Filtering criteria
perfect_F1 <- res$bayes.type == "Hybrid" &
              res$total >= 30 &
              res$max.p >= 0.9

F1_cells <- res[perfect_F1, ]

# Save perfect F1 barcodes
write.table(F1_cells, file = paste0(name, "_perfect_F1_cells.txt"), sep="\t", quote=F, row.names=FALSE)

###################################################################################################
## Step 5: Scatter plots (only total reads >= 30)
###################################################################################################
# Subset for plotting
res_plot <- res[res$total >= 30, ]

res_plot$bayes.type <- factor(res_plot$bayes.type, levels = c("P1", "Hybrid", "P2"))

base_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(linewidth = 0.5, color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5)
  )

max_value <- max(max(res_plot$P1_reads), max(res_plot$P2_reads))
res_plot_nona <- res_plot[!is.na(res_plot$bayes.type), ]

p1 <- ggplot(res_plot_nona, aes(x = P1_reads, y = P2_reads, color = bayes.type)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("P1" = "#1874CD", "Hybrid" = "#BFBFBF", "P2" = "#CD2626"),
    labels = c("P1-like", "Hybrid", "P2-like")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "P1 reads", y = "P2 reads", color = NULL,
       title = "Reads Scatter Plot (total reads ≥ 30)") +
  coord_fixed(xlim = c(0, max_value), ylim = c(0, max_value)) +
  base_theme

pdf(paste0(name, "_reads_scatter.pdf"), width=5.5, height=5)
print(p1)
dev.off()
