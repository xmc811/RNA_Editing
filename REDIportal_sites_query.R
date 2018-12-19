
setwd("/Users/mingchuxu/Documents/bioinformatics/RNA-edting")

library(tibble)
library(reshape2)

db_raw <- read.table("db_table", sep = "\t", stringsAsFactors = F, header = T)
decoy_raw <- read.table("decoy_table", sep = "\t", stringsAsFactors = F, header = T)


db_edit <- dcast(db_raw, Region.Position ~ Sample)
rm(db_raw)

decoy_edit <- dcast(decoy_raw, Region.Position ~ Sample)
rm(decoy_raw)

# read data from rna-editing analysis
# matrix ready for analysis


db_edit_frac <- c()

for (i in 1:nrow(db_edit)) {
        
        p <- 1 - sum(is.na(as.numeric(db_edit[i,2:ncol(db_edit)])))/(ncol(db_edit)-1)
        
        db_edit_frac <- append(db_edit_frac, p)
        
}

decoy_edit_frac <- c()

for (i in 1:nrow(decoy_edit)) {
        
        p <- 1 - sum(is.na(as.numeric(decoy_edit[i,2:ncol(decoy_edit)])))/(ncol(decoy_edit)-1)
        
        decoy_edit_frac <- append(decoy_edit_frac, p)
        
}


db_edit_med <- c()

for (i in 1:nrow(db_edit)) {
        
        p <- median(as.numeric(db_edit[i,2:ncol(db_edit)]), na.rm = T)
        
        db_edit_med <- append(db_edit_med, p)
        
}

decoy_edit_med <- c()

for (i in 1:nrow(decoy_edit)) {
        
        p <- median(as.numeric(decoy_edit[i,2:ncol(decoy_edit)]), na.rm = T)
        
        decoy_edit_med <- append(decoy_edit_med, p)
        
}



f_tier <- seq(from = 0.05, to = 0.8, by = 0.01)

f_decoy_stat <- c()

for (i in 1:length(f_tier)) {
        
        f_decoy_stat <- append(f_decoy_stat, sum(decoy_edit_frac>f_tier[i]))
        
}

ratio <- 1 - f_decoy_stat/f_stat

fraction_data <- data.frame(Fraction = f_tier, Datbase = f_stat, Decoy = f_decoy_stat, Ratio=ratio)

#fraction_data <- melt(fraction_data, id="Fraction")

library(ggplot2)

ggplot(data=fraction_data, aes(x=Fraction)) + 
        geom_line(aes(y=Real, color="Database")) +
        geom_line(aes(y=Decoy, color="Decoy")) +
        geom_line(aes(y=Ratio*1500, color="Ratio")) +
        scale_y_continuous(sec.axis = sec_axis(~./15, name = "Ratio [%]")) +
        ggtitle("Exonic A-to-I Editing Fractions in TCGA samples") + 
        geom_text(aes(label=ifelse(Fraction %in% c(0.05,0.1,0.2,0.3,0.4,0.5,0.60,0.70,0.8),Real,""),y=Real),hjust=-0.2, vjust=0) +
        geom_text(aes(label=ifelse(Fraction %in% c(0.05,0.1,0.2,0.3,0.4,0.5,0.60,0.70,0.8),Decoy,""),y=Decoy),hjust=-0.2, vjust=0) +
        scale_colour_brewer(palette="Set1")





