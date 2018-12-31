
RP_RE <- REDI_sites[,c("Chr", "Pos", "name", "strand", "name2", "endPos")]

get_pos_name <- function(chr, pos) {
        
        posName <- paste(chr, pos, sep = "-")
        
        return(posName)
}

for (i in 1:nrow(RP_RE)) {
        
        RP_RE$posName[i] <- get_pos_name(RP_RE$Chr[i], RP_RE$Pos[i])
        
}

get_window <- function(str, pos) {
        
        if (is.na(pos) == T) {
                return(NA)
        }
        
        Len5 <- pos[[1]][1]
        Len3 <- pos[[1]][4]
        
        if (str == "+") {
                return(list(c(Len5, Len3)))
        }  else {
                return(list(c(Len3, Len5)))
        }
        
}

for (i in 1:nrow(RP_RE)) {
        
        RP_RE$window[i] <- get_window(RP_RE$strand[i], RP_RE$endPos[i])
        
}

site_num <- function(ind) {
        
        return(sum(RP_RE$name2 == RP_RE$name2[ind]))

}

for (i in 1:nrow(RP_RE)) {
        
        RP_RE$siteNum[i] <- site_num(i)
        
}

ana <- RP_RE

edit_sample <- "SRR1562544"
edit_file <- "final_pair_table_brain"


edata <- read.table(edit_file, header = T, stringsAsFactors = F, sep = '\t')
edata <- subset(edata, edata[[1]] == edit_sample)
        
for (i in 1:nrow(ana)) {
                
        if (ana$posName[i] %in% edata[[2]]) {
                        
                ind <- match(ana$posName[i], edata[[2]])
                        
                ana$Edit[i] <- edata$Frequency[ind]
                        
        } else {
                        
                ana$Edit[i] <- NA
        }
                
} 

rp_file <- "Ribo_SRR1562539.Data"

data <- read.table(rp_file, header = F, stringsAsFactors = F, sep = '\t')
        
for (i in 1:nrow(ana)) {
                
        if (is.na(ana$endPos[i]) == T ) {
                        
                ana$Prof[i] <- NA
                        
        } else {
                        
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / sum * len
                        
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                        
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                        
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                        
                ana$Prof[i] <- list(sig)
                        
        }
}
        



