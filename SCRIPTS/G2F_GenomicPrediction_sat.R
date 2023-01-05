
library("dplyr")
library("agricolae")
library("ggplot2")
library("orthopolynom")
library("reshape2")
library("readxl")
library("lme4")
library("AGHmatrix", verbose=FALSE)
library('stringr')


pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############################################[Genomic Prediction: Reduced Model]############################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####GEBVs
data_reduced <- data.frame()

for(f in 1:12){
dirname <- "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/5.GenomicPrediction/Solutions"
filelist <- list.files(dirname)
filename <- file.path(dirname, filelist[f])
file.exists(filename)
filename
date = substr(filename, 112, 115)
UK = substr(filename, 116, 125)

solutions <- read.table(filename, header=T)  
head(solutions)  
summary(solutions)

solutions1 <- subset(solutions, solutions$effect==2 | solutions$effect==3 | solutions$effect== 4)
sol <- solutions1 %>%
  arrange(level)
head(sol)

pheno_dirname <- "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/5.GenomicPrediction/GenPredPheno"
pheno_filelist <- list.files(pheno_dirname)
pheno_filename <- file.path(pheno_dirname, pheno_filelist[f])
pheno_filename
file.exists(pheno_filename)
pheno <- read.table(pheno_filename, header=T)
T_mat <-pheno %>%
  dplyr::select(LocYear, EnvGradStandardized, leg0, leg1, leg2)
head(T_mat)
T_mat1 <- T_mat %>%
  distinct(LocYear, .keep_all=T) %>%
  arrange(., EnvGradStandardized)
head(T_mat1)

T_mat1 <- as.matrix(T_mat1[,c(3:5)])


LocYearEnvGradStandardized <- pheno %>%
  dplyr::select(LocYear, EnvGradStandardized) %>%
  distinct(LocYear, .keep_all=T) %>%
  arrange(., EnvGradStandardized)

for(i in 1:2126){
  sol2 <- subset(sol, sol$level==i)
  sol3 <- as.vector(sol2[,4])
  
  GEBV <- T_mat1%*%sol3
  
  Geno <- rep(i, 86)
  
  data6 <- cbind(Geno, f, EnvCov, GEBV)
  data_reduced <- rbind(data_reduced, data6)
  }
}

colnames(data_reduced) <- c("Geno", "SOLFILE", "LocYear", "EnvGradStandardized", "GEBV")
head(data_reduced)
summary(data_reduced)
genoID <-  read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/2.ModelSelection/geno2.txt_XrefID", header=F) 
colnames(genoID) <- c("Geno", "PedShort")

data_reduced <- merge(genoID, data_reduced, by="Geno") #Now I can match to phenotype.
head(data_reduced,20)
2126*12*86 #1462688. Worked for all 2126 hybrids, 8 GP cross validation(ish), 86 environments

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###############################################[Correlation/Bias BLUP]#####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
GenPredCor <- data.frame()
LocYearList <- unique(pheno$LocYear)

pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)

for(f in 1:length(unique(data_reduced$SOLFILE))){
   data10 <- subset(data_reduced, data_reduced$SOLFILE==f)

for(i in 1:length(unique(pheno$LocYear))){
  LY <- LocYearList[i]
  PhenoPerLocYear <- subset(pheno, pheno$LocYear==LY)
    head(PhenoPerLocYear)
  m1 <- lmer(GrainYield_Mg_ha_SpatialCor ~ (1|PedShort) + as.factor(Replicate), data=PhenoPerLocYear)
  
  blup = coef(m1)$PedShort
  blup$PedShort = rownames(blup)
  blup$BLUP = blup$`(Intercept)`
  blup = blup[,c(ncol(blup)-1,ncol(blup))]

  
  GEBVPerLocYear <- subset(data10, data10$LocYear==LY)


  GenPred <- merge(blup, GEBVPerLocYear, by="PedShort")    
  
  #Pearson's correlation
  cor_test <- cor(GenPred$GEBV, GenPred$BLUP)
  n = dim(GenPred)[1]

  vr = (1 - cor_test^2)/(n-2)
  rw = (cor_test/vr)/(1/vr)

  #Bias
  m2 <- lm(BLUP ~ GEBV, data=GenPred)
  plot(GenPred$GEBV, GenPred$BLUP)
  b1 = m2$coefficients[2]
  
  Intermediate <- cbind(f, GenPred[1,5], GenPred[1,6], cor_test, n, rw, b1)
  colnames(Intermediate) <- c("SOLFILE", "LocYear", "EnvGradStandardized", "Correlation", "Obs", "Correlation_Adjusted", "Bias")
  
  GenPredCor <- rbind(GenPredCor, Intermediate)
  }
}

head(GenPredCor)
86*12
for(i in 1:1032){
  if(GenPredCor[i,1] == 1){
    GenPredCor[i,8] = 2014
  } else if(GenPredCor[i,1] == 2){
    GenPredCor[i,8] = 2014
  } else if(GenPredCor[i,1] == 3){
    GenPredCor[i,8] = 2014
  } else if(GenPredCor[i,1] == 4){
    GenPredCor[i,8] = 2015
  } else if(GenPredCor[i,1] == 5){
    GenPredCor[i,8] = 2015
  } else if(GenPredCor[i,1] == 6){
    GenPredCor[i,8] = 2015
  } else if(GenPredCor[i,1] == 7){
    GenPredCor[i,8] = 2016
  } else if(GenPredCor[i,1] == 8){
    GenPredCor[i,8] = 2016
  } else if(GenPredCor[i,1] == 9){
    GenPredCor[i,8] = 2016
  } else if(GenPredCor[i,1] == 10){
    GenPredCor[i,8] = 2017
  } else if(GenPredCor[i,1] == 11){
    GenPredCor[i,8] = 2017
  } else if(GenPredCor[i,1] == 12){
    GenPredCor[i,8] = 2017
  }
}

for(i in 1:1032){
  if(GenPredCor[i,1] == 1){
    GenPredCor[i,9] = "Known"
  } else if(GenPredCor[i,1] == 2){
    GenPredCor[i,9] = "Unknown"
  } else if(GenPredCor[i,1] == 3){
    GenPredCor[i,9] = "NoGeno"
  } else if(GenPredCor[i,1] == 4){
    GenPredCor[i,9] = "Known"
  } else if(GenPredCor[i,1] == 5){
    GenPredCor[i,9] = "Unknown"
  } else if(GenPredCor[i,1] == 6){
    GenPredCor[i,9] = "NoGeno"
  } else if(GenPredCor[i,1] == 7){
    GenPredCor[i,9] = "Known"
  } else if(GenPredCor[i,1] == 8){
    GenPredCor[i,9] = "Unknown"
  } else if(GenPredCor[i,1] == 9){
    GenPredCor[i,9] = "NoGeno"
  } else if(GenPredCor[i,1] == 10){
    GenPredCor[i,9] = "Known"
  } else if(GenPredCor[i,1] == 11){
    GenPredCor[i,9] = "Unknown"
  } else if(GenPredCor[i,1] == 12){
    GenPredCor[i,9] = "NoGeno"
  }
}

colnames(GenPredCor) <- c("SOLFILE", "LocYear", "EnvGradStandardized", "Correlation", "Obs", "Correlation_Adjusted", "Bias", "Year.x","Type")

GenPredCor[c("Location", "Year")] <- str_split_fixed(GenPredCor$LocYear, "_", 2)

for(i in 1:1032){
  if(GenPredCor[i,8] != GenPredCor[i,11]){
    GenPredCor[i,4] =  NA
    GenPredCor[i,6] = NA
  }
}



plot(GenPredCor$Correlation, GenPredCor$Correlation_Adjusted)
GenPredCor$Correlation <- as.numeric(GenPredCor$Correlation)
GenPredCor$Bias <- as.numeric(GenPredCor$Bias)
GenPredCor <- na.omit(GenPredCor)
summary(GenPredCor)
head(GenPredCor)






GenPredCor_Known <- subset(GenPredCor, GenPredCor$Type=="Known")
GenPredCor_Known <- merge(GenPredCor_Known, data4, by="LocYear")
cor.test(GenPredCor_Known$herit, GenPredCor_Known$Correlation)
plot(GenPredCor_Known$herit, GenPredCor_Known$Correlation)

GenPredCor_Unknown <- subset(GenPredCor, GenPredCor$Type=="Unknown")
GenPredCor_Unknown <- merge(GenPredCor_Unknown, data4, by="LocYear")
cor.test(GenPredCor_Unknown$herit, GenPredCor_Unknown$Correlation)
plot(GenPredCor_Unknown$herit, GenPredCor_Unknown$Correlation)

GenPredCor_NoGeno <- subset(GenPredCor, GenPredCor$Type=="NoGeno")
GenPredCor_NoGeno <- merge(GenPredCor_NoGeno, data4, by="LocYear")
cor.test(GenPredCor_NoGeno$herit, GenPredCor_NoGeno$Correlation)
plot(GenPredCor_NoGeno$herit, GenPredCor_NoGeno$Correlation)


GenPredCor_BLUP <- GenPredCor
data_summary <- GenPredCor_BLUP %>%
  group_by(Type) %>%
  summarise(Mean_Cor = mean(Correlation, na.rm=T), 
            SD_Cor = sd(Correlation, na.rm=T),
            Min_cor = min(Correlation, na.rm=T),
            Max_cor = max(Correlation, na.rm=T),
            Mean_Bias = mean(Bias, na.rm=T),
            SD_Bias = sd(Bias, na.rm=T),
            Min_Bias = min(Bias, na.rm=T),
            Max_Bias = max(Bias, na.rm=T))
        
  

data_summary
data_summary <- GenPredCor_BLUP %>%
  group_by(SOLFILE) %>%
  summarise(Mean_Cor = mean(Correlation, na.rm=T), 
            SD_Cor = sd(Correlation, na.rm=T),
            Mean_Bias = mean(Bias, na.rm=T),
            SD_Bias = sd(Bias, na.rm=T))

data_summary$SOLKNOWN <- rep(c(1,2,3), 4)
data_summary$Year <- rep(c(2014, 2015, 2016, 2017), each=3)
data_summary

#write.csv(data_summary, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/5.GenomicPrediction/GP_GEBV.BLUP_sat.csv", quote = F, row.names = F)

ggplot(data_summary, aes(x=Year, y=Mean_Cor, fill=as.factor(SOLKNOWN))) +
  geom_bar(stat='identity', width=.85, position=position_dodge(width=0.9), colour="black", alpha=.9) +
  geom_errorbar(aes(ymin=Mean_Cor-SD_Cor, ymax=Mean_Cor+SD_Cor), width=.2, position=position_dodge(width=0.9), size=1.1) +
  ylim(0,1.04) +
  theme_bw() + 
  scale_fill_manual(labels = c("1" = "Known",
                               "2" = "Predicted",
                               "3" = "New Hybrids"),
                    values =c("cornflowerblue", "royalblue4", "midnightblue")) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  #theme(plot.margin = margin(.75, .75, .75, .75, "cm")) +
  theme(legend.box.background = element_rect(colour = "black")) +
  theme(legend.position = "none") +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Year",
       y = "Prediction Accuracy (r)",
       fill="Environmental Gradient")

plot(GenPredCor$EnvGradStandardized, GenPredCor$Correlation)


GenPredCorKnown = subset(GenPredCor, GenPredCor$SOLFILE==1 | GenPredCor$SOLFILE==3 | GenPredCor$SOLFILE==5 | GenPredCor$SOLFILE==7)

pheno_merge <- pheno[,c(2,8)]
pheno_merge <- unique(pheno_merge)
GenPredCorKnown1 <- merge(GenPredCorKnown, pheno_merge, by="LocYear")

summary(GenPredCorKnown1)
ggplot(data=GenPredCorKnown1, aes(x=reorder(LocYear,EnvGradStandardized.y,mean,na.rm=TRUE), y=Correlation)) +
  geom_point(col="royalblue4", na.rm=TRUE) +
  geom_smooth(col="red", method="gam")+
  ylim(-0.3, 0.7)+
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y=expression(bold(paste("Grain yield (kg h",a^{-1},")"))))


GenPredCorPredicted = subset(GenPredCor, GenPredCor$SOLFILE==2 | GenPredCor$SOLFILE==4 | GenPredCor$SOLFILE==6 | GenPredCor$SOLFILE==8)

GenPredCorPredicted1 <- merge(GenPredCorPredicted, pheno_merge, by="LocYear")
summary(GenPredCorPredicted1)
ggplot(data=GenPredCorPredicted1, aes(x=reorder(LocYear,EnvGradStandardized.y,mean,na.rm=TRUE), y=Correlation)) +
  geom_point(col="royalblue4", na.rm=TRUE) +
  geom_point(data=col="cornflowerblue", na.rm=TRUE) +
  ylim(-0.3, 0.7)+
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y=expression(bold(paste("Grain yield (kg h",a^{-1},")"))))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###############################################[Correlation Full Model]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

###Full Model BLUPF90 Solutions
solutions <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/3.GEBVs/solutions_leg2exp", header=T)  
head(solutions)  
summary(solutions)

solutions1 <- subset(solutions, solutions$effect==2 | solutions$effect==3 | solutions$effect== 4)
sol <- solutions1 %>%
  arrange(level)
head(sol)

pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)
T_mat <-pheno[,c(2,8:11)]
head(T_mat)
T_mat1 <- T_mat %>%
  distinct(LocYear, .keep_all=T) %>%
  arrange(., EnvGradStandardized)
head(T_mat1)

T_mat1 <- as.matrix(T_mat1[,c(3:5)])
LocYearEnvGradStandardized <- pheno[,c(2,8)]
EnvCov <- LocYearEnvGradStandardized %>%
  distinct(LocYear, .keep_all=T) %>%
  arrange(., EnvGradStandardized)

data_full <- data.frame()

for(i in 1:2126){
  sol2 <- subset(sol, sol$level==i)
  sol3 <- as.vector(sol2[,4])
  
  GEBV_full <- T_mat1%*%sol3
  Geno <- rep(i, 86)
  
  data6 <- cbind(Geno, EnvCov, GEBV_full)
  data_full <- rbind(data_full, data6)
}


LocYearList <- unique(pheno$LocYear)

GenPredCor <- data.frame()

for(f in 1:length(unique(data_reduced$SOLFILE))){
  data_reduced1 <- subset(data_reduced, data_reduced$SOLFILE==f)

data_cor <- data_full %>% 
  right_join(data_reduced1, by=c("LocYear","Geno"))

data_cor[c("Location", "Year")] <- str_split_fixed(data_cor$LocYear, "_", 2)

for(i in 1:length(unique(data_cor$LocYear))){

  LY <- unique(LocYearList)[i]
  GenPred <- subset(data_cor, data_cor$LocYear==LY)
  
  #Subset to just Hybrids that were actually testing in the environment
  pheno_int <- subset(pheno, pheno$LocYear==LY)
  PedShort <- as.data.frame(unique(pheno_int$PedShort))
  colnames(PedShort) <- c("PedShort")
  GenPred <- merge(PedShort, GenPred, by="PedShort")
  
  #Pearson's correlation
  cor_test <- cor(GenPred$GEBV, GenPred$GEBV_full)
  #Bias
  m2 <- lm(GEBV_full ~ GEBV, data=GenPred)
  plot(GenPred$GEBV, GenPred$GEBV_full)
  b1 = m2$coefficients[2]

  Intermediate <- cbind(f, GenPred[1,3], GenPred[1,4], cor_test, b1)
  colnames(Intermediate) <- c("SOLFILE", "LocYear", "EnvGradStandardized", "Correlation", "Bias")
  
  GenPredCor <- rbind(GenPredCor, Intermediate)
  }
}

head(GenPredCor)
86*12
for(i in 1:1032){
  if(GenPredCor[i,1] == 1){
    GenPredCor[i,6] = 2014
  } else if(GenPredCor[i,1] == 2){
    GenPredCor[i,6] = 2014
  } else if(GenPredCor[i,1] == 3){
    GenPredCor[i,6] = 2014
  } else if(GenPredCor[i,1] == 4){
    GenPredCor[i,6] = 2015
  } else if(GenPredCor[i,1] == 5){
    GenPredCor[i,6] = 2015
  } else if(GenPredCor[i,1] == 6){
    GenPredCor[i,6] = 2015
  } else if(GenPredCor[i,1] == 7){
    GenPredCor[i,6] = 2016
  } else if(GenPredCor[i,1] == 8){
    GenPredCor[i,6] = 2016
  } else if(GenPredCor[i,1] == 9){
    GenPredCor[i,6] = 2016
  } else if(GenPredCor[i,1] == 10){
    GenPredCor[i,6] = 2017
  } else if(GenPredCor[i,1] == 11){
    GenPredCor[i,6] = 2017
  } else if(GenPredCor[i,1] == 12){
    GenPredCor[i,6] = 2017
  }
}

for(i in 1:1032){
  if(GenPredCor[i,1] == 1){
    GenPredCor[i,7] = "Known"
  } else if(GenPredCor[i,1] == 2){
    GenPredCor[i,7] = "Unknown"
  } else if(GenPredCor[i,1] == 3){
    GenPredCor[i,7] = "NoGeno"
  } else if(GenPredCor[i,1] == 4){
    GenPredCor[i,7] = "Known"
  } else if(GenPredCor[i,1] == 5){
    GenPredCor[i,7] = "Unknown"
  } else if(GenPredCor[i,1] == 6){
    GenPredCor[i,7] = "NoGeno"
  } else if(GenPredCor[i,1] == 7){
    GenPredCor[i,7] = "Known"
  } else if(GenPredCor[i,1] == 8){
    GenPredCor[i,7] = "Unknown"
  } else if(GenPredCor[i,1] == 9){
    GenPredCor[i,7] = "NoGeno"
  } else if(GenPredCor[i,1] == 10){
    GenPredCor[i,7] = "Known"
  } else if(GenPredCor[i,1] == 11){
    GenPredCor[i,7] = "Unknown"
  } else if(GenPredCor[i,1] == 12){
    GenPredCor[i,7] = "NoGeno"
  }
}
head(GenPredCor)
colnames(GenPredCor) <- c("SOLFILE", "LocYear", "EnvGradStandardized", "Correlation", "Bias", "Year.x","Type")

GenPredCor[c("Location", "Year")] <- str_split_fixed(GenPredCor$LocYear, "_", 2)

for(i in 1:1032){
  if(GenPredCor[i,6] != GenPredCor[i,9]){
    GenPredCor[i,4] =  NA
  }
}

GenPredCor$SOLFILE <- as.character(GenPredCor$SOLFILE)
GenPredCor$Correlation <- as.numeric(GenPredCor$Correlation)
GenPredCor$Bias <- as.numeric(GenPredCor$Bias)
GenPredCor <- na.omit(GenPredCor)
summary(GenPredCor)

head(GenPredCor)
86*3

GenPredCor_Known <- subset(GenPredCor, GenPredCor$Type=="Known")
GenPredCor_Known <- merge(GenPredCor_Known, data4, by="LocYear")
cor.test(GenPredCor_Known$herit, GenPredCor_Known$Correlation)
plot(GenPredCor_Known$herit, GenPredCor_Known$Correlation)

GenPredCor_Unknown <- subset(GenPredCor, GenPredCor$Type=="Unknown")
GenPredCor_Unknown <- merge(GenPredCor_Unknown, data4, by="LocYear")
cor.test(GenPredCor_Unknown$herit, GenPredCor_Unknown$Correlation)
plot(GenPredCor_Unknown$herit, GenPredCor_Unknown$Correlation)

GenPredCor_NoGeno <- subset(GenPredCor, GenPredCor$Type=="NoGeno")
GenPredCor_NoGeno <- merge(GenPredCor_NoGeno, data4, by="LocYear")
cor.test(GenPredCor_NoGeno$herit, GenPredCor_NoGeno$Correlation)
plot(GenPredCor_NoGeno$herit, GenPredCor_NoGeno$Correlation)




GenPredCor_FULL <- GenPredCor
data_summary <- GenPredCor_FULL %>%
  group_by(Type) %>%
  summarise(Mean_Cor = mean(Correlation, na.rm=T), 
            SD_Cor = sd(Correlation, na.rm=T),
            Min_cor = min(Correlation, na.rm=T),
            Max_cor = max(Correlation, na.rm=T),
            Mean_Bias = mean(Bias, na.rm=T),
            SD_Bias = sd(Bias, na.rm=T),
            Min_Bias = min(Bias, na.rm=T),
            Max_Bias = max(Bias, na.rm=T))

data_summary


data_summary <- GenPredCor_FULL %>%
  group_by(SOLFILE) %>%
  summarise(Mean_Cor = mean(Correlation, na.rm=T), 
            SD_Cor = sd(Correlation, na.rm=T),
            Mean_Bias = mean(Bias, na.rm=T),
            SD_Bias = sd(Bias, na.rm=T))

data_summary$SOLKNOWN <- rep(c(1,2,3), 4)
data_summary$Year <- rep(c(2014, 2015, 2016, 2017), each=3)
data_summary
#write.csv(data_summary, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/5.GenomicPrediction/GP_GEBV.GEBVFULL_OnlyGenosInField_sat.csv", quote = F, row.names = F)

ggplot(data_summary, aes(x=Year, y=Mean_Cor, fill=as.factor(SOLKNOWN))) +
  geom_bar(stat='identity', width=.85, position=position_dodge(width=0.9), colour="black", alpha=.9) +
  geom_errorbar(aes(ymin=Mean_Cor-SD_Cor, ymax=Mean_Cor+SD_Cor), width=.2, position=position_dodge(width=0.9), size=1.1) +
  ylim(0,1.15) +
  theme_bw() + 
  scale_fill_manual(labels = c("1" = "Known",
                               "2" = "Predicted",
                               "3" = "New Hybrids"),
                    values =c("royalblue4", "cornflowerblue", "midnightblue")) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  #theme(plot.margin = margin(.75, .75, .75, .75, "cm")) +
  theme(legend.box.background = element_rect(colour = "black")) +
  theme(legend.position = "none") +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Year",
       y = "Prediction Accuracy (r)",
       fill="Environmental Gradient")



