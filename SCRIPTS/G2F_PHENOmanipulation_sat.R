#######This is a script to standardize the G2F yield data from 2014 - 2017

######Author: Seth Tolley
######Date: 20220623


######Needed Libraries
library("dplyr")
library("ggplot2")
library("orthopolynom")
library("reshape2")
library("usmap")
library("caret")
library("kernlab")
library("rnoaa")
library("soilDB")
library("stringr")
library("readxl")

####Import data
G2F14 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2014_hybrid_data_clean.csv", header=T)
G2F15 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2015_hybrid_data_clean.csv", header=T)
G2F16 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2016_hybrid_data_clean.csv", header=T)
G2F17 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2017_hybrid_data_clean.csv", header=T)


G2F <- rbind(G2F14, G2F15, G2F16, G2F17)

####Combining data across years
G2F$LocYear <- paste0(G2F$Field.Location,"_",G2F$Year)
G2F$PedLocYear <- paste0(G2F$Pedigree, "_",G2F$Field.Location,"_",G2F$Year)
length(unique(G2F$LocYear))
length(unique(G2F$Pedigree))

summary(G2F)
length(unique(G2F$LocYear))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#########################[How many LocYears were each hybrid tested in?]###########################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

###Final decision based on analysis in this area was to remove
###hybrids evaluated in less than 10 environments.
uniquePedLocYear <- G2F %>%
  dplyr::select(Pedigree, PedLocYear) %>%
  distinct(Pedigree, PedLocYear)
head(uniquePedLocYear)

uniqped <- unique(uniquePedLocYear$Pedigree)
data2 <- data.frame()
for(i in 1:length(uniqped)){
  
  SUM <- sum(uniquePedLocYear$Pedigree==uniqped[i])
  name <- uniqped[i]
  
  data1 <- cbind(name, SUM)
  data2 <- rbind(data2, data1)
  
}
data2$SUM <- as.numeric(data2$SUM)

summary(data2$SUM)
ggplot(data=data2, aes(x=SUM))+ 
  geom_histogram(binwidth=1, fill="steelblue", col="black") + 
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text = element_text(size=14)) +
  labs(x="Environments Each Hybrid Tested",
       y="Count")





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##########################[Correct Yield for Plot Size and Moisture]###############################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
G2F <- G2F[,c(39,40,1,2,5:10,15,31:34)]
colnames(G2F) <- c("LocYear", "PedLocYear", "Year", "Location","Pedigree",
                      "Replicate", "Block", "Plot", "Range", "Pass", "PlotArea",
                      "GrainMoisture", "TestWeight_lb.bu", "PlotWeight_lb", "GrainYield_bu.ac")

G2F <- G2F %>%
  filter(GrainYield_bu.ac>0)
length(unique(G2F$LocYear))


G2F$GrainYield_Mg_ha <- G2F$GrainYield_bu.ac*0.06277 #Convert bu/acre to Mg/ha
172*0.06277
14200000000*.0254

ggplot(data=G2F, aes(x=reorder(LocYear,GrainYield_Mg_ha,mean,na.rm=TRUE), y=GrainYield_Mg_ha)) +
  geom_boxplot(col="black", fill="steelblue4") +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y="Grain yield (kg/ha)")

uniqLocYear <- unique(G2F$LocYear)
data2 <- data.frame()

for(i in 1:length(uniqLocYear)){
  LY <- uniqLocYear[i]
  data1 <- subset(G2F, G2F$LocYear==LY)
  LocYearAvgYield <- mean(data1$GrainYield_kg.ha_15p, na.rm=TRUE)
  LocYear <- LY
  data1 <- cbind(LocYear, LocYearAvgYield)
  data2 <- rbind(data2, data1)
}

data2$LocYearAvgYield <- as.numeric(data2$LocYearAvgYield)




length(unique(G2F$LocYear)) #103 LocYears
ggplot(data=G2F, aes(x=reorder(LocYear,GrainYield_Mg_ha,mean,na.rm=TRUE), y=GrainYield_Mg_ha)) +
  geom_boxplot(col="black", fill="dodgerblue4", na.rm=TRUE) +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y=expression(bold(paste("Grain Yield (Mg h",a^{-1},")"))))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#################################[Estimate Repeatability]##############################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
list_LocYear <- unique(G2F$LocYear)

drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Repeatability_data <- data.frame()


for(j in 1:length(list_LocYear)){
  LE = list_LocYear[j]
  int1 <- subset(G2F, G2F$LocYear==LE)
  
  int1 <- int1 %>%
    drop_na(GrainYield_Mg_ha)
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  if(length(unique(int1$Replicate==1))){
    model <- lmer(GrainYield_Mg_ha ~ (1|Pedigree), data=int1) 
  } else {
    model <- lmer(GrainYield_Mg_ha ~ (1|Pedigree) +as.factor(Replicate), data=int1) 
  }
  summary(int1)
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) 
  
  RepNum <- length(unique(int1$Replicate))
  
  #deleting columns in variance component dataframe 
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  varComp
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  Repeatability <- sigmaG/(sigmaG + (sigmaE/length(unique(int1$Replicate))))
  
  #add columns to existing dataframe
  Data3 <- cbind(LE, RepNum, sigmaG, sigmaE, Repeatability)
  Repeatability_data <- rbind(Repeatability_data, Data3)
}
colnames(Repeatability_data) <- c("LocYear", "RepNum", "sigmaG", "sigmaE", "Repeatability")
Repeatability_noSpatial <- Repeatability_data
Repeatability_noSpatial %>%
  arrange(Repeatability)

###Perhaps we should remove the LocYears with a Repeatability less than .1 as an arbitrary cut off?
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#########################[Remove Plots Outside 3 SD of Mean]###########################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
list_LocYear <- unique(G2F$LocYear)

Outliers_data <- data.frame()

for(j in 1:length(list_LocYear)){
  LE = list_LocYear[j]
  int1 <- subset(G2F, G2F$LocYear==LE)
  int1 <- int1 %>%
    drop_na(GrainYield_Mg_ha)
  
  M = mean(int1$GrainYield_Mg_ha, na.rm=TRUE)
  SD = sd(int1$GrainYield_Mg_ha, na.rm=TRUE)
  
  for(i in 1:dim(int1)[1]){
    if(int1[i,16] > M+(3*SD)){
      int1[i,17] = 1
    } else if(int1[i,16] < M-(3*SD)){
      int1[i,17] = 1
    } else{
      int1[i,17] = 0
    }
  }
  Outliers_data <- rbind(Outliers_data, int1)
  print(j)
} 



Outliers_data$OutlierLogical <- Outliers_data$V17
summary(Outliers_data)
OnlyOutliers <- subset(Outliers_data, Outliers_data$OutlierLogical==1)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################[Field Spatial Correction]#############################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
OutliersRM_data <- subset(Outliers_data, Outliers_data$OutlierLogical==0)
OutliersRM_data <- OutliersRM_data %>%
  drop_na(GrainYield_Mg_ha)
length(unique(OutliersRM_data$LocYear))
list_LocYear <- unique(OutliersRM_data$LocYear)

Data3 <- data.frame()
SpATSRepeatability_data <- data.frame()
int4 <- data.frame()

pdf(file="~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/PLOTS/SpATS_SpatialCorrection.pdf",
    width=8,
    height=8)

for(j in 1:length(list_LocYear)){
  LE = list_LocYear[j]
  int1 <- subset(G2F, G2F$LocYear==LE)
  int2 <- int1 %>%
    drop_na(GrainYield_Mg_ha, Range, Pass)
  
  if(dim(int2)[1]>0){
    int2$R <- as.factor(int2$Pass)
    int2$C <- as.factor(int2$Range)
    
    m0 <- SpATS(response="GrainYield_Mg_ha", spatial = ~SAP(Range, Pass, nseg=c(10,20), degree=3, pord=2),
                genotype="Pedigree", genotype.as.random=T, random = ~ R + C, data=int2, 
                control= list(tolerance=1e-03))
    
    SpATS_Repeatability <- getHeritability(m0)
    
    Data3 <- cbind(LE, SpATS_Repeatability)
    SpATSRepeatability_data <- rbind(SpATSRepeatability_data, Data3)
    
    
    par(mar = c(1, 1, 1, 1)) # Set the margin on all sides to 2
    plot(m0)
    mtext(LE, 3, cex=1.2)
    
    int2$GrainYield_Mg_ha_SpatialCor <- predict(m0, newdata=int2)$predicted.values
    
    int3 <- int2 %>%
      dplyr::select(LocYear,Year,Replicate, Plot, Pedigree,GrainYield_Mg_ha_SpatialCor)
  } else {
    int3 <- int1 %>%
      dplyr::select(LocYear, Year, Replicate, Plot, Pedigree, GrainYield_Mg_ha)
    colnames(int3) <- c("LocYear", "Year", "Replicate", "Plot", "Pedigree", "GrainYield_Mg_ha_SpatialCor")
  }
  int4 <- rbind(int4, int3)
  print(j) #Status
}
dev.off()


SpATSRepeatability_data %>%
  arrange(SpATS_Repeatability)

SpATS_Correction <- int4

myCol_BLUE <- c("darkslategray1", "lightsteelblue2", "royalblue1", "midnightblue")
pdf(file="~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/PLOTS/G2F_YieldPerLocYear_W.WOOutliers.pdf",
    width=30,
    height=10)

#No outlier removal
ggplot(data=G2F, aes(x=reorder(LocYear,GrainYield_Mg_ha,mean,na.rm=TRUE), y=GrainYield_Mg_ha, fill=as.factor(Year))) +
  geom_boxplot(col="black", na.rm=TRUE) +
  scale_fill_manual(values = myCol_BLUE) + 
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y=expression(bold(paste("Grain yield (Mg h",a^{-1},")"))),
       fill="Year",
       title="All Phenotypic Data")

#Outliers in red 
ggplot(data=Outliers_data, aes(x=reorder(LocYear,GrainYield_Mg_ha,mean,na.rm=TRUE), y=GrainYield_Mg_ha, 
                               fill=as.factor(Year))) +
  geom_boxplot(col="black", na.rm=TRUE) +
  scale_fill_manual(values = myCol_BLUE) + 
  geom_point(data=OnlyOutliers, aes(x=reorder(LocYear,GrainYield_Mg_ha,mean,na.rm=TRUE), y=GrainYield_Mg_ha),col="red", pch=4, size=3)+
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  guides(col="none") +
  labs(x="Location Years",
       y=expression(bold(paste("Grain yield (Mg h",a^{-1},")"))),
       fill="Year",
       title="Outliers (More than 3 SD from Mean) X'd Out ")


ggplot(data=SpATS_Correction, aes(x=reorder(LocYear,GrainYield_Mg_ha_SpatialCor,mean,na.rm=TRUE), y=GrainYield_Mg_ha_SpatialCor, 
                                  fill=as.factor(Year))) +
  geom_boxplot(col="black", na.rm=TRUE) +
  scale_fill_manual(values = myCol_BLUE) + 
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  guides(col="none") +
  labs(x="Location Years",
       y=expression(bold(paste("Grain yield (Mg h",a^{-1},")"))),
       fill="Year")
       #title="SpATS Spatial Correction After Outlier Removal")

dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#################[Estimate Repeatability After Spatial Correction]#####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
drops <- c("var1","var2","sdcor") 
Data3 <- data.frame()
Repeatability_data <- data.frame()
head(SpATS_Correction)

for(j in 1:length(list_LocYear)){
  LE = list_LocYear[j]
  int1 <- subset(SpATS_Correction, SpATS_Correction$LocYear==LE)
  
  int1 <- int1 %>%
    drop_na(GrainYield_Mg_ha_SpatialCor)
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  if(length(unique(int1$Replicate.x==1))){
    model <- lmer(GrainYield_Mg_ha_SpatialCor ~ (1|Pedigree), data=int1) 
  } else {
    model <- lmer(GrainYield_Mg_ha_SpatialCor ~ (1|Pedigree) +as.factor(Replicate), data=int1) 
  }
  summary(int1)
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) 
  
  RepNum <- length(unique(int1$Replicate))

  
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  varComp
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  Repeatability <- sigmaG/(sigmaG + (sigmaE/length(unique(int1$Replicate))))
  
  #add columns to existing dataframe
  Data3 <- cbind(LE, RepNum, sigmaG, sigmaE, Repeatability)
  Repeatability_data <- rbind(Repeatability_data, Data3)
}
colnames(Repeatability_data) <- c("LocYear", "RepNum_SC", "sigmaG_SC", "sigmaE_SC", "Repeatability_SC")
Repeatability_data <- Repeatability_data %>%
  arrange(Repeatability)


Repeatability_data <- merge(Repeatability_noSpatial, Repeatability_data, by="LocYear")
plot(Repeatability_data$Repeatability, Repeatability_data$Repeatability_SC)
abline(a=0, b=1)


LocYeartoRemove = unique(c("NYH1_2015", "TXH2_2015","NEH1_2016", "NEH4_2016", "NYH1_2016", "SCH1_2016", "TXH2_2016", "ILH1_2017", "ILH2_2017", "INH1_2017", "TXH2_2017", #No lat/long
                    "MNH1_2017", "NYH2_2016", "NYH3_2016", "COH1_2017", #Bad phenotypic data from note in metadata
                    "ARH1_2016", "ARH2_2016", "GAH1_2016", #Irrigated but no amount given
                    "ARH1_2017", "ILH1_2017", "NEH4_2017", "ARH2_2017", "NEH1_2014" #Repeatability less than 0.1
                    ))

data1  <- subset(SpATS_Correction, !(SpATS_Correction$LocYear %in% LocYeartoRemove)) #Removing the 17 Environments with a repeatability less than 0.1.

genoswithG2F <- unique(data1$Pedigree)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###################################[Where were the hybrids tested?]################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)


ggplot(data=pheno, aes(x=reorder(LocYear,Year,mean,na.rm=TRUE), y=PedShort, fill=GrainYield_kg.ha_15p)) +
  geom_tile(fill="black") +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  labs(x="Location Years",
       y="Hybrid")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
########################[Standardizing Environmental Gradient from -1 to 1]########################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####I also made the second order legendre polynomial for RRM

data2 <- data1 %>%
  drop_na(GrainYield_Mg_ha_SpatialCor)
model1 <- lmer(GrainYield_Mg_ha_SpatialCor ~ 0 + (1|Pedigree) + LocYear, data=data1)


EG <- t(data.frame(as.list(fixef(model1))))
EG <- cbind(rownames(EG), data.frame(EG, row.names=NULL))



colnames(EG) <- c("L1", "EnvGrad")
EG$L1 <- gsub("LocYear", "", as.character(EG$L1))
colnames(EG) <- c("LocYear", "EnvGrad")

EG[EG == "TXH1.Early_2017"] <- "TXH1-Early_2017"; EG[EG == "TXH1.Late_2017"] <- "TXH1-Late_2017"; EG[EG == "TXH1.Dry_2017"] <- "TXH1-Dry_2017"

data2 <- merge(data1, EG, by="LocYear")
length(unique(data2$LocYear))



G2F.min <- min(data2$EnvGrad,na.rm=TRUE)
G2F.max <- max(data2$EnvGrad,na.rm=TRUE)

data2$EnvGradStandardized <- -1+2*((data2$EnvGrad-G2F.min)/(G2F.max-G2F.min))

leg4coef <- legendre.polynomials(n=2, normalized=T)
leg4 <- as.matrix(as.data.frame(polynomial.values(polynomials=leg4coef, x=data2$EnvGradStandardized)))

colnames(leg4) <- c("leg0", "leg1", "leg2")
leg4 <- leg4[, 1:ncol(leg4)]
dim(leg4)              
                  
data2 <- cbind(data2, leg4)
                

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#######################################[Genotypic Data]############################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

setwd("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/GenotypicInfo/")

G2F_Geno <- read.table("G2F_Imputed_2173HybridGeno_100878.txt", header=T, sep="\t")
G2F_Geno[1:5,1:5]

geno1 <- G2F_Geno[,1]
length(unique(geno1))

geno2 <- G2F_Geno[,-1]
geno2[geno2==1] <- 2
geno2[geno2==0.5] <- 1

G2F_Geno <- cbind(geno1, geno2)
G2F_Geno[1:5, 1:5]
colnames(G2F_Geno)[1] <- "Taxa"
length(unique(G2F_Geno$Taxa))

#This dataset can be used to get the Taxa into the phenotypic Pedigree column
GenoToPheno <- read_excel("g2f_2017_gbs_hybrid_codes copy.xlsx", sheet=2)
head(GenoToPheno)
GenoToPheno$PedShort <- str_pad(GenoToPheno$PedShort, 6, pad="0")
head(GenoToPheno)

G2F_Geno1 <- merge(GenoToPheno, G2F_Geno, by="Taxa")
G2F_Geno1[1:5,1:5]

data3 <- merge(GenoToPheno, data2, by="Pedigree")
head(data3)

PedShort <- as.data.frame(unique(data3$PedShort))
PedShort$PedShort <- PedShort$`unique(data3$PedShort)`
PedShort <- PedShort[,-1]
PedShort <- as.data.frame(PedShort)
G2F_Geno2 <- merge(PedShort, G2F_Geno1, by="PedShort")
G2F_Geno3 <- G2F_Geno2[,-c(2:3)]
G2F_Geno3[1:5,1:5]
G2F_Geno3[1:5,100870:100879]


write.table(G2F_Geno3, file="G2F_geno.txt", quote=F, sep = " ", row.names = F)

G2F_Geno3[1:5,1:5]




PedShort <- as.data.frame(unique(G2F_Geno3$PedShort))
PedShort$PedShort <- PedShort$`unique(G2F_Geno3$PedShort)`
PedShort <- PedShort[,-1]
PedShort <- as.data.frame(PedShort)
G2F_pg <- merge(PedShort, data3, by="PedShort")
G2F_pg <- G2F_pg[,-c(2:3)]
head(G2F_pg)
write.table(G2F_pg, file="G2F_pheno1.txt", quote=F, sep = " ", row.names = F)
length(unique(G2F_pg$PedShort))


save.image(file="~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/SCRIPTS/G2F_PHENOmanipulationWorkspace_sat.RData")




G2F_env <- G2F_pg %>%
  group_by(LocYear) %>%
  summarize(meanYield = mean(GrainYield_Mg_ha_SpatialCor, na.rm=TRUE))

G2F_year <- G2F_pg %>%
  group_by(Year) %>%
  summarize(meanYield = mean(GrainYield_Mg_ha_SpatialCor, na.rm=TRUE))
G2F_year







G2F <- G2F %>%
  filter(GrainYield_bu.ac>0)

G2F_summary <- G2F %>%
  filter(LocYear %in% LocYeartoRemove)
length(unique(G2F_summary$LocYear))

G2F_summary <- subset(G2F, !(LocYear %in% c(LocYeartoRemove)))


ONH2_2017 <- subset(G2F_summary, G2F_summary$LocYear=="ONH2_2017")
notONH2_2017<- subset(G2F_summary, G2F_summary$LocYear!="ONH2_2017")

summary(ONH2_2017) #38 days to anthesis
summary(notONH2_2017) #70
