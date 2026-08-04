###################################################################################################
## Define function to plot marker accessibility scores 
###################################################################################################
plot.act.scores    <- function(df, 
                               acts=acts, 
                               info=NULL, 
                               top=NULL,
                               logT=F,
                               marker.dist=NULL,
                               outname="markerActivityScores.png", 
                               lim=0.95){

    # prep data
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,which(rownames(df) %in% colnames(acts))]

    # reorder rows
    rownames(info) <- info$geneID
    info <- info[order(info$type),]
    info.genes <- rownames(info)
    act.genes <- rownames(acts)
    rd.cells <- rownames(df)

    # common genes
    common <- intersect(info.genes, act.genes)
    info <- info[which(rownames(info) %in% common),]
    info.ordered <- rownames(info)
    sub.scores <- acts[info.ordered, , drop = FALSE]
    gids <- info.ordered

    # -------------------- # Use a 1x1 layout with a smaller canvas when plotting a single gene; otherwise use the multi-gene layout --------------------
    if (length(gids) == 1) {
        png(file=outname, width=4, height=4, units="in", res=500, type="cairo")
        layout(matrix(1, ncol=1, nrow=1))
    } else {
        nrows <- ceiling(length(gids)/8)
        totals <- nrows*8
        ratio <- nrows/8

        png(file=outname, width=12, height=ratio*12, units="in", res=500, type="cairo") 
        layout(matrix(c(1:totals), ncol=8, byrow=T))
    }
    # -----------------------------------------------------------------------

    par(mar=c(2,2,1,1))

    # adjust cluster IDs
    message("begin plotting pre-defined markers...")
    for (i in 1:length(gids)){

        # copy meta data
        gene.index <- match(gids[i], rownames(sub.scores))
        acv <- as.numeric(sub.scores[gene.index, , drop = TRUE])

        # set up plot cols/sizes
        orderRow <- order(acv, decreasing=F)
        #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
        #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
        #cols <- inferno(100)
        #cols <- plasma(100)
        cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.75)(100)
        acv <- as.numeric(acv[orderRow])
        if(logT==T){
            acv <- log2(acv+1)
        }
        df2 <- df[orderRow,]
        acv[is.na(acv)] <- 0
        acv[is.infinite(acv)] <- 0
        upper.lim <- quantile(acv, lim)
        acv[acv > upper.lim] <- upper.lim
        if(!is.null(marker.dist)){
            message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
            colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
        }else{
            min.acv <- min(acv) - (1e-6*min(acv))
            max.acv <- max(acv) + (1e-6*max(acv))
            message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
	    if(min.acv == max.acv){
		next
	    }
            colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
        }
        colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
	colvec[is.na(colvec) & acv == 0] <- cols[1]
        sizes <- 0.2 #rescale(acv, c(0.25, 0.3))

        plot(df2$umap1, df2$umap2, col=colvec,
             main=info$name[i], # No gene name on the plot, or if wat to keep it, then set main=info$name[i]
             xlab="", ylab="", bty="n",
             xaxt="n", yaxt="n", pch=16, cex=0.2)

    }

    # turn device off
    dev.off()

}

#--------------------------------------------------------------------------------------
#                        plot selected marker accessibility scores
#---------------------------------------------------------------------------------------
library(Matrix)
library(RColorBrewer) 

# ================== Parse command-line arguments ==================
# Usage: Rscript plot_markers.R <rds_file> <marker_file> <gene_arg> <output_dir> <output_prefix>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_markers.R <rds_file> <marker_file> <gene_arg> <output_dir> <output_prefix>")
}

rds_file <- args[1]  
marker_file <- args[2]
gene_arg <- args[3]  # # e.g. "tb1,gt1,GRAS17" or a file
output_dir <- args[4]  # Directory to save plots
output_prefix<- args[5]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =================== load data ===========================
a <- readRDS(rds_file)
acts <- a$impute.activity   # plot imputed activity scores

df <- data.frame(
  umap1 = a$b$umap1,
  umap2 = a$b$umap2,
  row.names = rownames(a$b)
)

info <- read.table(marker_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (file.exists(gene_arg)) {
  # case: input is a file with one gene name per line
  names_vec <- scan(gene_arg, what = "character", sep = "\n", quiet = TRUE)
} else {
  # case: input is a comma-separated string
  names_vec <- trimws(strsplit(gene_arg, ",")[[1]])
}


# Keep only cells that are present in both 'df' and 'acts'
common.cells <- intersect(rownames(df), colnames(acts))
df <- df[common.cells, , drop = FALSE]
acts <- acts[, common.cells, drop = FALSE]

# ======================= Select markers by name ==========================
names_vec <- unique(names_vec)
info_sub  <- subset(info, name %in% names_vec & geneID %in% rownames(acts))
info_sub  <- info_sub[match(names_vec, info_sub$name), ] 
rownames(info_sub) <- info_sub$geneID

max_panels <- 180   # maximum number of genes per figure to avoid oversized PNG (Cairo error)
chunks <- split(seq_len(nrow(info_sub)), (seq_len(nrow(info_sub)) - 1) %/% max_panels)

# ========================================== Plot selected markers=========================================
part <- 1
for (idx in chunks) {
  info_part <- info_sub[idx, , drop = FALSE]
  out_file  <- file.path(output_dir, sprintf("%s_combined.part%02d.png", output_prefix, part))
  plot.act.scores(
    df     = df,
    acts   = acts,
    info   = info_part,
    logT   = FALSE,
    outname= out_file,
    lim    = 0.999
  )
  message("Saved: ", out_file)
  part <- part + 1
}

