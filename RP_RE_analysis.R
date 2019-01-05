library(ggplot2)


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

edit_sample_b1 <- "SRR1562544"
edit_sample_b2 <- "SRR1562545"
edit_sample_b3 <- "SRR1562546"
edit_sample_b4 <- "SRR1562547"
edit_sample_b5 <- "SRR1562548"

edit_file <- "final_pair_table_brain"


edata <- read.table(edit_file, header = T, stringsAsFactors = F, sep = '\t')


edata_b1 <- subset(edata, edata[[1]] == edit_sample_b1)
edata_b2 <- subset(edata, edata[[1]] == edit_sample_b2)
edata_b3 <- subset(edata, edata[[1]] == edit_sample_b3)
edata_b4 <- subset(edata, edata[[1]] == edit_sample_b4)
edata_b5 <- subset(edata, edata[[1]] == edit_sample_b5)
        
for (i in 1:nrow(ana)) {
                
        if (ana$posName[i] %in% edata_b1[[2]]) {
                        
                ind <- match(ana$posName[i], edata_b1[[2]])
                        
                ana$EditB1[i] <- edata_b1$Frequency[ind]
                        
        } else {
                        
                ana$EditB1[i] <- NA
        }
                
} 

for (i in 1:nrow(ana)) {
        
        if (ana$posName[i] %in% edata_b2[[2]]) {
                
                ind <- match(ana$posName[i], edata_b2[[2]])
                
                ana$EditB2[i] <- edata_b2$Frequency[ind]
                
        } else {
                
                ana$EditB2[i] <- NA
        }
        
} 

for (i in 1:nrow(ana)) {
        
        if (ana$posName[i] %in% edata_b3[[2]]) {
                
                ind <- match(ana$posName[i], edata_b3[[2]])
                
                ana$EditB3[i] <- edata_b3$Frequency[ind]
                
        } else {
                
                ana$EditB3[i] <- NA
        }
        
} 

for (i in 1:nrow(ana)) {
        
        if (ana$posName[i] %in% edata_b4[[2]]) {
                
                ind <- match(ana$posName[i], edata_b4[[2]])
                
                ana$EditB4[i] <- edata_b4$Frequency[ind]
                
        } else {
                
                ana$EditB4[i] <- NA
        }
        
} 

for (i in 1:nrow(ana)) {
        
        if (ana$posName[i] %in% edata_b5[[2]]) {
                
                ind <- match(ana$posName[i], edata_b5[[2]])
                
                ana$EditB5[i] <- edata_b5$Frequency[ind]
                
        } else {
                
                ana$EditB5[i] <- NA
        }
        
} 

rm(edata_b1)
rm(edata_b2)
rm(edata_b3)
rm(edata_b4)
rm(edata_b5)







rp_file_b1 <- "Ribo_SRR1562539.Data"
rp_file_b2 <- "Ribo_SRR1562540.Data"
rp_file_b3 <- "Ribo_SRR1562541.Data"
rp_file_b4 <- "Ribo_SRR1562542.Data"
rp_file_b5 <- "Ribo_SRR1562543.Data"


data <- read.table(rp_file_b1, header = F, stringsAsFactors = F, sep = '\t')
        
for (i in 1:nrow(ana)) {
                
        if (is.na(ana$endPos[i]) == T ) {
                        
                ana$ProfB1[i] <- NA
                        
        } else {
                        
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / (sum+1) * len
                        
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                        
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                        
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                        
                ana$ProfB1[i] <- list(sig)
                        
        }
}


data <- read.table(rp_file_b2, header = F, stringsAsFactors = F, sep = '\t')

for (i in 1:nrow(ana)) {
        
        if (is.na(ana$endPos[i]) == T ) {
                
                ana$ProfB2[i] <- NA
                
        } else {
                
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / (sum+1) * len
                
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                
                ana$ProfB2[i] <- list(sig)
                
        }
}


data <- read.table(rp_file_b3, header = F, stringsAsFactors = F, sep = '\t')

for (i in 1:nrow(ana)) {
        
        if (is.na(ana$endPos[i]) == T ) {
                
                ana$ProfB3[i] <- NA
                
        } else {
                
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / (sum+1) * len
                
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                
                ana$ProfB3[i] <- list(sig)
                
        }
}


data <- read.table(rp_file_b4, header = F, stringsAsFactors = F, sep = '\t')

for (i in 1:nrow(ana)) {
        
        if (is.na(ana$endPos[i]) == T ) {
                
                ana$ProfB4[i] <- NA
                
        } else {
                
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / (sum+1) * len
                
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                
                ana$ProfB4[i] <- list(sig)
                
        }
}


data <- read.table(rp_file_b5, header = F, stringsAsFactors = F, sep = '\t')

for (i in 1:nrow(ana)) {
        
        if (is.na(ana$endPos[i]) == T ) {
                
                ana$ProfB5[i] <- NA
                
        } else {
                
                site_data <- subset(data, data$V4 == ana$posName[i])
                sum <- sum(site_data$V6)
                len <- ana$window[i][[1]][1] + ana$window[i][[1]][2]
                sig <- site_data$V6 / (sum+1) * len
                
                if (ana$strand[i] == "-") {
                        sig <- rev(sig)
                }
                
                LenUp <- ana$window[i][[1]][1]
                LenDown <- ana$window[i][[1]][2]
                
                sig <- c(rep(NA,bin_around-LenUp),sig,rep(NA,bin_around-LenDown))
                
                ana$ProfB5[i] <- list(sig)
                
        }
}






plot_data <- data.frame(cbind(seq(from = -500, to = 500, by = 1), subset(ana, posName=="chr19-20807628")[[11]][[1]]))

ggplot(plot_data, aes(x=X1, y=X2)) + geom_line() + 
        geom_area(fill = "steelblue1") +
        geom_vline(xintercept = 0, colour = "firebrick")

a <- c()
b <- c()

summary_data <- subset(ana, EditB5<0.2 & EditB5<1.4 & is.na(ProfB5)==F & siteNum<=400)

summary_ribo <- function(data) {
        
        df <- data.frame()
        
        for (i in 1:nrow(data)) {
                
                 df <- rbind(df, data$ProfB5[i][[1]])
        }
        
        vec <- c()
        
        for (i in 1:1001) {
              
                val <- mean(df[,i], na.rm = T)
                vec <- append(vec, val)
        }
        
        plot_df <- data.frame(cbind(seq(from = -500, to = 500, by = 1), vec))
        
        colnames(plot_df) <- c("X1", "X2")
        
        ggplot(plot_df, aes(x=X1, y=X2)) + geom_line(size = 2) +
                geom_area(fill = "steelblue1") +
                geom_vline(xintercept = 0, colour = "firebrick") +
                geom_hline(yintercept = 1, colour = "black")
        return(vec)
}

b <- append(b, summary_ribo(summary_data))

a <- append(a, summary_ribo(summary_data))



summary_site <- function(site) {
        
        data <- subset(ana, posName==site)
        
        df <- data.frame()
        
        for (i in 1:5) {
                
                df <- rbind(df, data[1,i+14][[1]])
        }
        
        vec <- c()
        
        for (i in 1:1001) {
                
                val <- mean(df[,i], na.rm = T)
                vec <- append(vec, val)
        }
        
        plot_df <- data.frame(cbind(seq(from = -500, to = 500, by = 1), vec))
        
        colnames(plot_df) <- c("X1", "X2")
        
        ggplot(plot_df, aes(x=X1, y=X2)) + geom_line(size = 1.5) +
                geom_area(fill = "dodgerblue4") +
                geom_vline(xintercept = 0, colour = "firebrick") +
                geom_hline(yintercept = 1, colour = "black") +
                #scale_y_continuous(limits = c(0,8)) +
                xlab("Position") +
                ylab("Ribosomal Occupancy")
        
}

summary_site("chr15-65425334")




#final

a_mat <- data.frame(matrix(a, nrow = 5, byrow = T))

b_mat <- data.frame(matrix(b, nrow = 5, byrow = T))

nm <- c()

for (i in 1:1001) {
        
        nm[i] <- mean(a_mat[,i])
        
}


ct <- c()

for (i in 1:1001) {
        
        ct[i] <- mean(b_mat[,i])
        
}

plot_data <- data.frame(cbind(seq(from = -500, to = 500, by = 1), nm))
plot_data_ct <- data.frame(cbind(seq(from = -500, to = 500, by = 1), ct))

ggplot(plot_data, aes(x=V1, y=ct)) + geom_line(size=1.5) + 
        geom_area(fill = "skyblue1") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_hline(yintercept = 1, colour = "black") +
        xlab("Position") +
        ylab("Ribosomal Occupancy") +
        scale_y_continuous(limits = c(0,2.5))
