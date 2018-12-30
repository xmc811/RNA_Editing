
setwd("/Users/mingchuxu/Documents/repos/RNA_Editing")

library(tibble)
library(dplyr)

RefSeq <- read.table("RefSeq_hg19.Data", header = T, sep = '\t', stringsAsFactors = F, na.strings = "")

RefSeq$TxLen <- with(RefSeq, abs(txEnd - txStart))
RefSeq$cdsLen <- with(RefSeq, abs(cdsEnd - cdsStart))

RefSeq <- RefSeq[with(RefSeq, order(name2, -cdsLen, -exonCount)),]

RefSeq <- RefSeq[!duplicated(RefSeq$name2),]

RefSeq <- subset(RefSeq, cdsLen > 0)

RefSeq <- RefSeq[with(RefSeq, order(chrom, txStart)),]


REDI_sites <- read.table("original_REDI_sites", header = F, sep = ':', stringsAsFactors = F)

colnames(REDI_sites) <- c("Chr", "Pos")

site_in_tx <- function(chr, pos) {
        
        data <- subset(RefSeq, chrom == chr & txStart <= pos & txEnd >= pos)
        
        nrow(data)
}

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$inTx[i] <- site_in_tx(REDI_sites$Chr[i], REDI_sites$Pos[i])
        
}

REDI_sites <- subset(REDI_sites, inTx > 0)



get_annotation <- function(chr, pos) {
        
        data <- subset(RefSeq, chrom == chr & txStart <= pos & txEnd >= pos)
        
        data[1,]
}

site_annotation <- data.frame()

for (i in 1:nrow(REDI_sites)) {
        
        site_annotation <- rbind(site_annotation, get_annotation(REDI_sites$Chr[i], REDI_sites$Pos[i]))
        
}

REDI_sites <- cbind(REDI_sites, site_annotation)
rm(site_annotation)

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$exonStart[i] <- list(as.numeric(strsplit(REDI_sites[i,13],",")[[1]]))
        REDI_sites$exonEnd[i] <- list(as.numeric(strsplit(REDI_sites[i,14],",")[[1]]))
        
}

get_exon_len <- function(lst1, lst2) {
        
        lst1[[1]] - lst2[[1]] + 1
        
}

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$exonLen[i] <- list(get_exon_len(REDI_sites[i,23], REDI_sites[i,22]))
        
}

find_exon_num <- function(pos, lst1, lst2) {
        
        mat <- matrix(data = append(lst1[[1]],lst2[[1]]), ncol = 2)
        
        for (i in 1:nrow(mat)) {
                
                if (pos >= mat[i,1] & pos <=mat[i,2]) {
                        
                        return(list(append(i, mat[i,2] - pos)))
                }
                
        }
        
        return(NA)
        
}

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$exonPos[i] <- find_exon_num(REDI_sites$Pos[i],REDI_sites$exonStart[i],REDI_sites$exonEnd[i])
        
}

bin_around = 500

incre_vec <- function(vec) {
        
        new_vec <- 0
        
        for (i in 1:length(vec)) {
                
                new_vec <- append(new_vec, sum(vec[1:i]))
                
        }
        return(new_vec)
}


find_end <- function(exonPos, exonLen) {
        
        if (is.na(exonPos) == T) {
                return(NA)
        }
        
        Num <- exonPos[[1]][1]
        Loc <- exonPos[[1]][2]
        
        fullLen <- sum(exonLen[[1]])
        
        LocLen <- sum(exonLen[[1]][1:Num]) - Loc
        
        Len3 <- min(fullLen - LocLen, bin_around)
        Len5 <- min(LocLen - 1, bin_around)
        
        End3 <- LocLen + Len3
        End5 <- LocLen - Len5
        
        vec <- incre_vec(exonLen[[1]])
        
        
        for (i in 1:(length(vec)-1)) {
                
                if ((End5 > vec[i]) & (End5 <= vec[i+1])) {
                        
                        Ex5 <- i
                        Loc5 <- vec[i+1] - End5
                        
                }
        }
        
        for (i in 1:(length(vec)-1)) {
                
                if ((End3 > vec[i]) & (End3 <= vec[i+1])) {
                        
                        Ex3 <- i
                        Loc3 <- vec[i+1] - End3
                        
                }
        }
        
        return(list(c(Len5, Ex5, Loc5, Len3, Ex3, Loc3)))
        
}

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$endPos[i] <- find_end(REDI_sites$exonPos[i],REDI_sites$exonLen[i])
        
}

find_bed_bound <- function(exStart, exEnd, endPos) {
        
        if (is.na(endPos) == T) {
                return(NA)
        }

        mat <- matrix(data = append(exStart[[1]],exEnd[[1]]), ncol = 2)
        
        bound <- c()
        
        Ex5 <- endPos[[1]][2]
        Ex3 <- endPos[[1]][5]
        
        Loc5 <- endPos[[1]][3]
        Loc3 <- endPos[[1]][6]
        
        for (i in Ex5:Ex3) {
                
                bound <- append(bound, mat[i,1])
                bound <- append(bound, mat[i,2])
                
        }
                
        bound[1] <- bound[2] - Loc5
        bound[length(bound)] <- bound[length(bound)] - Loc3
        
        return(list(bound))
        
}

for (i in 1:nrow(REDI_sites)) {
        
        REDI_sites$bedBound[i] <- find_bed_bound(REDI_sites$exonStart[i], REDI_sites$exonEnd[i], REDI_sites$endPos[i])
        
}

generate_bed <- function(ind) {
        
        if (is.na(REDI_sites$bedBound[ind]) == T) {
                
                bed <- data.frame()
                
                return(bed)
                
        }
        
        chr <- REDI_sites$Chr[ind]
        pos <- REDI_sites$Pos[ind]
        
        posName <- paste(chr, pos, sep = "-")
        
        gene <- REDI_sites$name2[ind]
        
        mat <- matrix(data = REDI_sites$bedBound[ind][[1]], ncol = 2, byrow = T)
        
        vec <- c()
        
        for (i in 1:nrow(mat)) {
                
                range <- seq(from = mat[i,1], to = mat[i,2])
                vec <- append(vec, range)
                
        }
        
        len <- length(vec)
        
        bed <- data.frame("Chr" = rep(chr, times = len), "Start" = vec, "End" = vec, "PosN" = rep(posName,times = len), "Gene" = rep(gene, times = len))
        
        return(bed)
        
}

all_bed <- data.frame()
all_bed2 <- data.frame()


for (i in 1:3000) {
        
        all_bed <- bind_rows(all_bed, generate_bed(i))
}

for (i in 3001:nrow(REDI_sites)) {
        
        all_bed2 <- bind_rows(all_bed2, generate_bed(i))
}

bedfile <- bind_rows(all_bed, all_bed2)

rm(all_bed)
rm(all_bed2)

write.table(bedfile, file = "ribo_profiling_bed.Data", quote = F, sep = "\t", row.names = F, col.names = F)

