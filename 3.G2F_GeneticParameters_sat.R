
library("dplyr")
library("ggplot2")
library("orthopolynom")
library("reshape2")
library("readxl")
library("lme4")
library("AGHmatrix", verbose=FALSE)
library("MASS")
library("scales")


pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)


ggplot(data=pheno, aes(x=reorder(LocYear,GrainYield_Mg_ha_SpatialCor,mean,na.rm=TRUE), y=GrainYield_Mg_ha_SpatialCor, 
                                  fill=as.factor(Year))) +
  geom_boxplot(col="black", na.rm=TRUE) +
  scale_fill_manual(values = myCol_BLUE) + 
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x=element_text(angle=90, size=12)) +
  guides(col="none") +
  labs(x="Location Years",
       y=expression(bold(paste("Grain Yield (Mg h",a^{-1},")"))),
       fill="Year")





T_mat <-pheno[,c(8:11)]
head(T_mat)
T_mat1 <- T_mat %>%
  distinct(EnvGradStandardized, .keep_all=T) %>%
  arrange(., EnvGradStandardized)
head(T_mat1)
T_mat1 <- as.matrix(T_mat1[,c(2:4)])

G_data <- as.matrix(   c(  0.80836,      0.24529,     -0.14625,    
                           0.24529,      0.15604,     -0.55873E-02,
                           -0.14625,     -0.55873E-02,  0.68149E-01))

G_mat <- matrix(G_data, nrow=3, ncol=3, byrow=T)
head(G_mat)



x <- T_mat1%*%G_mat
dim(x)
head(x)
dim(T_mat1)
y <- x%*%t(T_mat1)
dim(y)
head(T_mat)
EnvCov <- T_mat %>%
  distinct(EnvGradStandardized, .keep_all=F) %>%
  arrange(., EnvGradStandardized)

residual <- exp(0.53385 + EnvCov$EnvGradStandardized*-0.47724E-01)
plot(EnvCov$EnvGradStandardized,residual, xlab = "Standardized Environmental Gradient", ylab="Residual")
residual

###Heritability at each environment
data4 <- data.frame()
for(i in 1:86){
  sigmaA <- y[i,i]
  hetres <- residual[i]
  herit <- sigmaA/(sigmaA+hetres)
  data3 <- cbind(i, herit)
  data4 <- rbind(data4, data3)
}
data4 <- cbind(data4, EnvCov)
plot(data4$EnvGradStandardized, data4$herit)
head(data4)

ggplot(data4, aes(x=EnvGradStandardized, y=herit))+
  geom_point(size=2, col="royalblue4")+
  ylim(c(0,.4)) +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text=element_text(size=16)) +
  theme(legend.position = "right") +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Standardized Environmental Gradient",
       y = expression(bold(paste("Heritability (", h^2,")"))))



head(pheno)
###Genetic Correlation
LocYearEnvGradStandardized <- pheno[,c(2,8)]
EnvCov <- LocYearEnvGradStandardized %>%
  distinct(EnvGradStandardized, .keep_all=T) %>%
  arrange(., EnvGradStandardized)

data4 <- merge(data4, EnvCov, by="EnvGradStandardized")
data5 <- data.frame()
for(i in 1:86){
  for(j in 1:86){
    cov1 <- y[i,j]
    var1 <- y[i,i]
    var2 <- y[j,j]
    EnvGradStandardized1 <- EnvCov[i,2]
    EnvGradStandardized2 <- EnvCov[j,2]
    LocYear1 <- EnvCov[i,1]
    LocYear2 <- EnvCov[j,1]
    gen_cor <- cov1/sqrt(var1*var2)
    data3 <- cbind(i,j,EnvGradStandardized1, EnvGradStandardized2, LocYear1, LocYear2, gen_cor)
    data5 <- rbind(data5,data3)
  }
}

head(data5)
str(data5)
summary(data5)

data5$i <- as.numeric(data5$i)
data5$j <- as.numeric(data5$j)
data5$gen_cor <- as.numeric(data5$gen_cor)
summary(data5)

cols = c("white", "lavender", "powderblue", "cornflowerblue","royalblue2", "midnightblue")
ggplot(data5, aes(x=reorder(LocYear1,i, mean,na.rm=TRUE), y=reorder(LocYear2,j, mean,na.rm=TRUE), fill=gen_cor))+
#ggplot(data5, aes(x=EnvGradStandardized1, y=EnvGradStandardized2))+
  geom_tile()+
  theme_bw() +
  #scale_fill_viridis(discrete=F, option="cividis", direction=-1)
  #scale_fill_gradient2(low="wheat", mid="goldenrod3", high="black", midpoint=.2,
   #                    limits=c(-0.2, 1)) +
  scale_fill_gradientn(colors=cols,
                       limits=c(0, 1)) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(angle=90, size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(legend.position = "right") +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Location Years",
       y = "Location Years",
       fill="Correlation")

  
  
####GEBVs
solutions <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/3.GEBVs/solutions_leg2exp", header=T)  
head(solutions)  
summary(solutions)

solutions1 <- subset(solutions, solutions$effect==2 | solutions$effect==3 | solutions$effect== 4)
sol <- solutions1 %>%
  arrange(level)
head(sol)


data7 <- data.frame()
for(i in 1:2126){
sol2 <- subset(sol, sol$level==i)
sol3 <- as.vector(sol2[,4])

GEBV <- T_mat1%*%sol3

Geno <- rep(i, 86)


data6 <- cbind(Geno, EnvCov, GEBV)
data7 <- rbind(data7, data6)
}
head(data7)


GEBVs_to_col <- subset(data7, data7$EnvGradStandardized==1)
GEBVs_to_col1 <- GEBVs_to_col %>%
  arrange(GEBV, keep.all=TRUE)
GEBVs_to_col2 <- as.vector(GEBVs_to_col1[c(1:10, 2117:2126),1])
head(GEBVs_to_col2,20)
GEBVs_to_col3 <- data7[data7$Geno %in% GEBVs_to_col2,]
summary(data7)
cols = c("319" = "cornflowerblue", "462"="cornflowerblue", "452"="cornflowerblue", "270"="cornflowerblue",  "461"="cornflowerblue",
         "278"="cornflowerblue", "282"= "cornflowerblue",  "316"="cornflowerblue", "290"="cornflowerblue", "315"="cornflowerblue",
         "1522"="royalblue4", "494"="royalblue4", "1513"="royalblue4", "1499"="royalblue4", "1494"="royalblue4",
         "1246"="royalblue4", "1422"="royalblue4", "1501"="royalblue4", "1488"="royalblue4", "1502"="royalblue4")
ggplot(data=data7, aes(x=EnvGradStandardized, y=GEBV, group=Geno)) +
  #geom_point() +
  geom_smooth(data= data7, size=.04, se=F, col="grey50", alpha=.25) +
  geom_smooth(data=GEBVs_to_col3, linewidth=1, alpha=.5,
              aes(group=as.factor(Geno), col=as.factor(Geno))) +
  scale_colour_manual(values=cols) +
  guides(col="none") +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text=element_text(size=16)) +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Standardized Environmental Gradient",
       y = "Genomic Estimated Breeding Values")
  







####SNP Solutions
snp_eff <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/4.SNPEffect/snp_sol", header=TRUE)



head(snp_eff)
summary(snp_eff)


data9 <- data.frame()
for(i in 1:max(snp_eff$snp)){
snp_eff2 <- subset(snp_eff, snp_eff$snp==i)
snp_eff2 <- snp_eff2[,6]
SNP_Effect <- T_mat1%*%snp_eff2

snp <- rep(i, 86)

data8 <- cbind(snp, EnvCov, SNP_Effect)
data9 <- rbind(data9, data8)
}
 





t1 <- rep(1:86, 100878)
data9.1 <- cbind(data9, t1)


data9.1$AbsSNP_Effect = abs(data9.1$SNP_Effect)

MapFile <- read.table("~/Documents/Purdue/Phd/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_MapFile.txt", sep="\t", header=T)
MapFile$snp <- 1:100878

data9.1 <- merge(data9.1, MapFile, by="snp")
summary(data9.1)
 
data11 <- data.frame()

for(i in 1:86){
data9.2 <- subset(data9.1, data9.1$t1==i)

data9.3 <- data9.2 %>%
  arrange(desc(AbsSNP_Effect), keep.all=T)

data10 <- data9.3[c(1:30),]

data10.1 <- data10 %>%
  arrange(POS, keep.all=T)

for(j in 1:20){
  if(data10.1[j+1,9] - data10.1[j,9] <1000){ #Removing SNPs that are in LD with one another
    data10.1[j+1,6] = NA
  }
}

data10.2 <- data10.1 %>%
  arrange(desc(AbsSNP_Effect), keep.all=T)
data10.2 <- na.omit(data10.2)

data10.3 <- data10.2[c(1:10),]
data11 <- rbind(data11, data10.3)
}

head(data11,30)



summary(data11)

# distinct
Table_SNPs <- data11 %>% 
  distinct(snp, .keep_all = T)

Table_SNPs
head(data9.1)

wantedsnps <- as.vector(unique(data11$snp))
data9.4 <- data9.1[data9.1$snp %in% wantedsnps,]
length(unique(data9.4$snp))
write.csv(data11, "~/Desktop/table1.csv")
LocYearEnvGradStandardized <- pheno[,c(2,8)]
EnvCov <- LocYearEnvGradStandardized %>%
  distinct(EnvGradStandardized, .keep_all=T) %>%
  arrange(., EnvGradStandardized)


rows = c(1:86)
times = length(unique(data11$snp))
EnvCov1 <- EnvCov[rep(rows, times),]

data9.5 <- cbind(data9.4, EnvCov1)
data9.6 <- data9.5[,c(1:5)]
summary(data9.6)

MapFile <- read.table("~/Documents/Purdue/Phd/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_MapFile.txt", sep="\t", header=T)
MapFile$snp <- 1:100878

data9.7 <- merge(data9.6, MapFile, by="snp")
summary(data9.7)
data9.7 <- subset(data9.7, data9.7$SNP_ID!="S2_236187946")

cols = c("darkred", "red2", "indianred1", "pink", "white", "lavender", "cornflowerblue","royalblue2", "midnightblue")

ggplot(data9.7, aes(x=reorder(LocYear,EnvGradStandardized,mean,na.rm=TRUE), y=SNP_ID, fill=SNP_Effect))+
  geom_tile() +
  theme_bw() +
  #scale_fill_gradient2(low="wheat", mid="goldenrod3", high="black", midpoint=.9,
   #                    limits=c(0, 1.55)) +
    scale_fill_gradientn(colors=cols,
                         limits=c(-.002, .002)) +
  #scale_fill_gradient2(low="firebrick", mid="grey90", high="steelblue", midpoint=0,
   #                   limits=c(-1.5, 1.55)) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(angle=90, size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(legend.position = "right") +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Location Years",
       y = "SNP Identification",
       fill="SNP Effect")

head(data9.7)
cols = c("darkred", "red2", "indianred1", "pink", "white", "lavender", "cornflowerblue","royalblue2", "midnightblue")

cols = c("bisque1", "pink", "cornsilk3", "goldenrod1", "goldenrod3", "darkorange", "darkorange3", "indianred1", "red2", "darkred", 
         "springgreen1", "springgreen3", "forestgreen", "darkgreen", "cyan", "deepskyblue", "cornflowerblue", "royalblue4", "mediumblue", 
         "midnightblue", "grey55", "black")

ggplot(data9.7, aes(x=EnvGradStandardized, y=SNP_Effect, color=SNP_ID))+
#geom_point()+
  geom_line(stat="summary", fun=mean, linewidth=1) +
  theme_bw() +
  scale_color_manual(values=cols) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(size=14)) +
  theme(axis.text.y=element_text(size=14)) +
  theme(legend.position = "right") +
  guides(color=guide_legend(ncol=1)) +
  #scale_x_discrete(position="top") +
  #scale_y_discrete(position="right") +
  labs(x="Standardized Environmental Gradient",
       y=expression(bold(paste("SNP Effect of Major Allele (Mg h",a^{-1},")"))),
       color="SNP Identification")




head(data9.7)
unique(data9.7$snp)


data9.8 <- data9.7 %>%
  distinct(SNP_ID, .keep_all = T)


data9.8
TOPSNPS <- unique(data9.8$SNP_ID)
TOPSNPS <- c("PedShort", TOPSNPS)
TOPSNPS <- TOPSNPS[-4]






pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)

EnvCov <- pheno %>%
  dplyr::select(LocYear, EnvGrad) %>%
  distinct()
LocYearList <- unique(pheno$LocYear)

head(pheno)
pheno_BLUPs <- data.frame()

  for(i in 1:length(unique(pheno$LocYear))){
    LY <- LocYearList[i]
    PhenoPerLocYear <- subset(pheno, pheno$LocYear==LY)
    head(PhenoPerLocYear)
    m1 <- lmer(GrainYield_Mg_ha_SpatialCor ~ (1|PedShort) + as.factor(Replicate), data=PhenoPerLocYear)
    
    blup = coef(m1)$PedShort
    blup$PedShort = rownames(blup)
    blup$BLUP = blup$`(Intercept)`
    blup = blup[,c(ncol(blup)-1,ncol(blup))]
    
    data10 <- cbind(LY, blup)
    pheno_BLUPs <- rbind(pheno_BLUPs, data10)
    
  }

    head(pheno_BLUPs)
colnames(pheno_BLUPs) <- c("LocYear", "PedShort", "GrainYield_Mg_ha_SpatialCor_BLUP")

geno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_geno.txt", header=T)
geno[1:5,1:5]
geno_sub <- geno %>%
  dplyr::select(TOPSNPS)
geno_sub[1:5,1:5]
pheno1 <- merge(pheno_BLUPs, geno_sub, by="PedShort")
pheno2 <- merge(EnvCov, pheno1, by="LocYear")
head(pheno2)

List_LocYear = unique(pheno2$LocYear)


Data14 <- data.frame()
Data15 <- data.frame()

for(i in 1:length(List_LocYear)){
  for(j in 5:25){
    Data13 <- data.frame()
    LY = List_LocYear[i]
    int1 <- subset(pheno2, pheno2$LocYear==LY)

    EG = int1[1,2]
    SNP_ID = colnames(int1)[j]
    colnames(int1)[j] = "SNP"
    int1[,j] = as.factor(int1[,j])
    
    m1 = lm(GrainYield_Mg_ha_SpatialCor_BLUP ~ SNP, int1)
    anova(m1)
    
    m2 <- LSD.test(m1, trt=c("SNP"), DFerror=anova(m1)[2,1], MSerror = anova(m1[2,3]), alpha=0.05, group=TRUE)
    m3 <- m2$groups
      m3$MajorAllele = row.names(m2$groups)
    m3 <- m3 %>%
      arrange(MajorAllele)
    m3
    Data13[1,1] <- anova(m1)[1,5]
    Data13[1,2] <- if(Data13[1,1]<0.001){
      "***"
    } else if (Data13[1,1]<0.01){
      "**"
    } else if (Data13[1,1]<0.05){
      "*"
    }else {
      "NS"
    }
    m2$means
    d0 = m2$means[1,1]
    d0sd = m2$means[1,2]
    r0 = m2$means[1,3]
    Tukey0 = m3[1,2]
    
    d1 = m2$means[2,1]
    d1sd = m2$means[2,2]
    r1 = m2$means[2,3]
    Tukey1 = m3[2,2]
    
    d2 = m2$means[3,1]
    d2sd = m2$means[3,2]
    r2 = m2$means[3,3]
    Tukey2 = m3[3,2]
    
    Data14 <- cbind(LY, EG, SNP_ID, Data13, d0, d0sd,Tukey0, d1, d1sd, Tukey1, d2, d2sd, Tukey2, r0, r1, r2)
    Data15 <- rbind(Data15, Data14)
  }
}
head(Data15)
colnames(Data15) = c("LocYear", "EnvGrad", "SNP_ID", "ANOVA_Pval", "ANOVA_sig", "Mean_Minor", "SD_Minor", "Tukey_Minor",
                     "Mean_Het", "SD_Het", "Tukey_Het", "Mean_Major", "SD_Major", "Tukey_Major", "Num_Minor", "Num_Het", "Num_Major")  
22*86
head(Data15)

for(i in 1:dim(Data15)[1]){
if(Data15[i,5]=="NS"){
  Data15[i,8]= "NA";
  Data15[i,11]= "NA";
  Data15[i,14]= "NA"
  }
}
  
write.csv(Data15, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/7.SNPEffectDif/G2F_SNPEffectInDifferentEnvironments_sat.csv", quote = F,
          row.names = F)




geno_sub1 <- geno_sub
  
geno_sub2 <- subset(geno_sub1, geno_sub1$S2_20959401==2 &
                      geno_sub1$S2_236187932==0 &
                      geno_sub1$S3_214902774==2 &
                      geno_sub1$S3_215373678==0 &
                      geno_sub1$S3_215380499==0 &
                      geno_sub1$S3_216647852==2 &
                      geno_sub1$S3_216794651==2 &
                      geno_sub1$S3_219184259==2 &
                      geno_sub1$S3_222732398==2 &
                      #geno_sub1$S4_33885271== 2 
                      geno_sub1$S4_37012538 == 2 
                      #geno_sub1$S4_156482117 == 2 &
                      #geno_sub1$S4_156994628 == 2 &
                      )










data12 <- data9.7 %>%
  dplyr::select(snp, SNP_ID, CHR, POS)%>%
  distinct(snp, .keep_all = TRUE)

data12
#save.image(file = "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/SCRIPTS/G2F_GeneticParameters.RData")





LD <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/LinkageDisequilibrium/plink.ld", header=T)

str(LD)
LD$Dist_bp <- LD$BP_B - LD$BP_A
summary(LD)




myCol_BLUE <- c("#B4B4B4FF", "#A2ACB4FF", "#90A4B4FF", "#7E9BB4FF", "#6C93B4FF", "#5A8BB4FF", "#4883B4FF", "#367BB4FF", "#2473B4FF", "#0966B4FF")





ggplot(data=LD, aes(x=Dist_bp, y=R2, colour=as.factor(CHR_A))) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE, col='black', size=1.1) +
  geom_hline(yintercept=0.2, col="red", linetype='dashed') +
  scale_colour_manual(values = myCol_BLUE) + 
  ylim(0,1) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(size=14)) +
  theme(axis.text.y=element_text(size=14)) +
  theme(legend.box.background = element_rect(colour = "black")) +
  theme(legend.position = c(0.8,0.7)) +
  labs(x="Physical Distance",
       y=expression(bold("Linkage Disequilibrium ("~R^2~")")),
       colour="Chromosome")
      


eigenval <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/PCA/G2F_2126Hybrids_PCA.eigenval", header=F)
head(eigenval)
plot(eigenval$V1)


pc <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/RESULTS/PCA/G2F_2126Hybrids_PCA.eigenvec", header=F)
head(pc)
plot(pc$V3, pc$V4)
