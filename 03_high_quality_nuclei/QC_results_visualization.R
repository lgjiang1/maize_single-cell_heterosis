library(MASS)
library(viridis)
library(dplyr)  

# Arguments
args <- commandArgs(T)
if(length(args) != 2){ stop("Rscript step3_QC_results_visualization.R [input] [name]") } # Two arguments required
#Rscript step3_QC_results_visualization.R LSC3_BxK.updated_metadata.txt LSC3_BxK

# Load arguments
input <- as.character(args[1])
name <- as.character(args[2])

# Load data
message(" - loading meta data for ", name)
meta <- read.table(input, header=T, sep="\t")

# Ensure data types are correct
meta$pTSS <- as.numeric(meta$pTSS)
meta$tss_z <- as.numeric(meta$tss_z)
meta$FRiP <- as.numeric(meta$FRiP)
meta$acr_z <- as.numeric(meta$acr_z)
meta$pOrg <- as.numeric(meta$pOrg)
meta$gREF <- as.numeric(meta$gREF)
meta$bREF <- as.numeric(meta$bREF)

pdf(paste(name, ".QC_FIGURES.pdf", sep=""), width=12, height=3)
layout(matrix(c(1:5), nrow=1))

###### **Step 1: Filtering based on total insertions (barcode rank plot)** #######
plotDist <- function(x, main="") {
    x <- x[order(x$total, decreasing=T),]  # Sort by total insertions
    rank <- log10(seq(1, nrow(x)))
    depth <- log10(x$total + 1)  # Compute log10(total)
    df <- data.frame(rank=rank, depth=depth)

    # Filter cells with depth >= 10^3
    df.n <- subset(df, df$depth > 3)
    cells <- nrow(df.n)
    
    if (cells == 0) {
        warning("No cells remaining after filtering!")
        return(NULL)
    }

    knee <- rank[cells]
    reads <- as.integer(10^(df.n$depth[nrow(df.n)]))

    # Plot barcode rank distribution
    plot(rank[(cells+1):length(rank)], depth[(cells+1):length(rank)],
         type="l", lwd=2, col="grey75", main=main,
         xlim=range(rank),
         ylim=range(depth),
         xlab="Barcode rank (log10)",
         ylab="Unique Tn5 insertions (log10)")
    lines(rank[1:cells], depth[1:cells], lwd=2, col="darkorchid4")
    grid()
    abline(v=knee, col="red", lty=2, lwd=1)
    abline(h=log10(reads), col="red", lty=2, lwd=1)
    text(2,1, labels=paste("# cells=", cells, ", # reads=", reads-1, sep=""))

    # Return only the rows corresponding to filtered cells
    return(x[x$total >= 1000, ])  
}
#============Step 1: Filtering based on total insertions=======================
message(" - filtering cells based on total insertions ...")
meta.filtered1 <- plotDist(meta, main=name)
message(" - number of cells after total filtering = ", nrow(meta.filtered1))



#============Step 2: Filtering cells based on pTSS > 0.2 and tss_z > -2==========
message(" - filtering cells based on pTSS and tss_z ...")

meta.filtered2 <- meta.filtered1 %>% filter(pTSS > 0.2, tss_z > -2)
message(" - number of cells after pTSS and tss_z filtering = ", nrow(meta.filtered2))

# **Plot TSS QC**
den_pTSS <- kde2d(log10(meta.filtered1$total), meta.filtered1$pTSS, n=300, h=c(0.2, 0.05),
                   lims=c(c(2.5, 6.5), c(0,1.0)))

image(den_pTSS, useRaster=T, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="Fraction reads in TSS",
      main=name)
grid(lty=1, lwd=0.5, col="grey90")
abline(h=0.2, col="red", lty=2, lwd=1)
legend("topright", legend=paste("# cells = ", nrow(meta.filtered2), sep=""), fill=NA, col=NA, border=NA)
box()



#============Step 3: Filtering cells based on FRiP > 0.2 and acr_z > -2==============
message(" - filtering cells based on FRiP and acr_z ...")

meta.filtered3 <- meta.filtered2 %>% filter(FRiP > 0.2, acr_z > -2)
message(" - number of cells after FRiP and acr_z filtering = ", nrow(meta.filtered3))

# **Plot FRiP QC**
den_FRiP <- kde2d(log10(meta.filtered2$total), meta.filtered2$FRiP, n=300, h=c(0.2, 0.05),
                   lims=c(c(2.5, 6.5), c(0,1.0)))

image(den_FRiP, useRaster=T, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="Fraction Tn5 insertions in ACRs",
      main=name)
grid(lty=1, lwd=0.5, col="grey90")
abline(h=0.2, col="red", lty=2, lwd=1)
legend("topright", legend=paste("# cells = ", nrow(meta.filtered3), sep=""), fill=NA, col=NA, border=NA)
box()

#=================Step 4: Filtering cells based on pOrg < 0.05========================
message(" - filtering cells based on pOrg < 0.05 ...")

meta.filtered4 <- meta.filtered3 %>% filter(pOrg < 0.05)
message(" - number of cells after pOrg filtering = ", nrow(meta.filtered4))

# **Plot pOrg QC**
den_pOrg <- kde2d(log10(meta.filtered3$total), meta.filtered3$pOrg, n=300, h=c(0.2, 0.05),
                   lims=c(c(2.5, 6.5), c(0,1.0)))

image(den_pOrg, useRaster=T, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="Fraction reads in Organelles",
      main=name)
grid(lty=1, lwd=0.5, col="grey90")
abline(h=0.05, col="red", lty=2, lwd=1)
legend("topright", legend=paste("# cells = ", nrow(meta.filtered4), sep=""), fill=NA, col=NA, border=NA)
box()

#=============Step 5: Computing cut-off and filtering cells with `call == 1`==============
message(" - computing cut-off based on call == 0 ...")

# Compute the median of `dif` for cells where `call == 0`
meta.filtered4$dif <- meta.filtered4$gREF - meta.filtered4$bREF
cut.off <- median(meta.filtered4$dif[meta.filtered4$call == 0], na.rm = TRUE)
message(" - cut-off value = ", cut.off)

# Filtering `call == 1`
message(" - filtering cells based on call == 1 ...")
meta.filtered5 <- meta.filtered4 %>% filter(call == 1)
message(" - number of cells after call == 1 filtering = ", nrow(meta.filtered5))

# **Plot gREF - bREF QC**
den_gREF <- kde2d(log10(meta.filtered4$total), meta.filtered4$dif,
                   n=300, h=c(0.2, 0.05),
                   lims=c(c(2.5, 6.5), c(min(meta.filtered4$dif, na.rm=T), max(meta.filtered4$dif, na.rm=T))))

image(den_gREF, useRaster=T, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="gREF - bREF",
      main=name)
grid(lty=1, lwd=0.5, col="grey90")
abline(h = ifelse(cut.off < 0, 0, cut.off), col = "red", lty = 2, lwd = 1)  # Mark cut-off threshold
legend("topright", legend=paste("# cells = ", nrow(meta.filtered5), sep=""), fill=NA, col=NA, border=NA)
box()

write.table(meta.filtered5, file=paste(name, ".final_filtered_metadata.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
dev.off()

