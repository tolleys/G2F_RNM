#######This is a script to prepare pheno and geno files for BLUPF90

######Author: Seth Tolley
######Date: 20220706


######Needed Libraries
library("dplyr")
library('tidyr')
library("ggplot2")
library("orthopolynom")
library("reshape2")
library("readxl")
library("lme4")
library("stringr")
library("orthopolynom")
library("reshape2")
library("usmap")
library("caret")
library("soilDB")
library("readxl")
library("randomForest")
library("nasapower")
library("doParallel")


pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##########################[Cluster Locations Based on Lat+Long and Weather]########################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cnames <- c("Location", "City", "Latitude", "Longitude")
w14 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/EnvironmentalInfo/g2f_2014_field_characteristics.csv")[,c(3,4,8,7)]
w15 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/EnvironmentalInfo/g2f_2015_field_metadata.csv")[,c(1,4,6,7)]
w16 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/EnvironmentalInfo/g2f_2016_field_metadata.csv", header=F)[-1,c(1,3,9,10)]
w17 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/EnvironmentalInfo/g2f_2017_field_metadata.csv", header=F)[-1,c(1,3,9,10)]
colnames(w14) <- cnames; colnames(w15) <- cnames; colnames(w16) <- cnames; colnames(w17) <- cnames
w14$Year = 2014; w15$Year = 2015; w16$Year = 2016; w17$Year = 2017
head(w14); head(w15); head(w16); head(w17)

g2f_loc <- rbind(w14, w15, w16, w17)
g2f_loc$LocYear <- paste0(g2f_loc$Location,"_",g2f_loc$Year)
head(g2f_loc)
g2f_loc <- g2f_loc[(g2f_loc$LocYear %in% unique(pheno$LocYear)),]

locYear <- unique(pheno$LocYear)[!(unique(pheno$LocYear) %in% g2f_loc$LocYear)] #There are no LocYear that we don't have  
#phenotypic data


g2f_loc <- g2f_loc %>%
  mutate_all(na_if,"")
g2f_loc1 <- na.omit(g2f_loc)
str(g2f_loc1)
g2f_loc1$Latitude <- as.numeric(g2f_loc1$Latitude)
g2f_loc1$Longitude <- as.numeric(g2f_loc1$Longitude)




d1 <- pheno[,c(2,7)]
d100 <- unique(d1)
g2f_loc2 <- merge(g2f_loc1, d100, by="LocYear")

head(g2f_loc2)
AllCounty <- map_data("county")
AllState <- map_data("state")
ggplot() + 
  geom_polygon(data=AllCounty, aes(x=long, y=lat, group=group),
               color="grey80", fill="grey80", size = .1) +
  geom_polygon(data=AllState, aes(x=long, y=lat, group=group),
               color="white", fill="grey80",  size = .5, alpha = .2) +
  geom_point(data=g2f_loc2, aes(x=Longitude, Latitude, col=EnvGrad, 
                                shape=as.factor(Year), size=as.factor(Year))) +
  scale_shape_manual(values=c(1, 2, 3, 4)) +
  scale_color_gradient2(low="firebrick",mid="cornflowerblue", high="darkblue", midpoint=6,
                        limits=c(0, 12)) +
  theme_void() +
  labs(shape="Year",
       color="Average Yield")



####Import data for planting and harvesting date.
G2F14 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2014_hybrid_data_clean.csv", header=T)
G2F15 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2015_hybrid_data_clean.csv", header=T)
G2F16 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2016_hybrid_data_clean.csv", header=T)
G2F17 <- read.csv("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/HybridPhenotypes/g2f_2017_hybrid_data_clean.csv", header=T)
  G2F <- rbind(G2F14, G2F15, G2F16, G2F17)
  G2F$LocYear <- paste0(G2F$Field.Location,"_",G2F$Year)
  G2F$PedLocYear <- paste0(G2F$Pedigree, "_",G2F$Field.Location,"_",G2F$Year)

g2f_plantdate <- G2F %>%
  dplyr::select(LocYear, Date.Planted, Date.Harvested) %>%
  distinct(LocYear, .keep_all = T)
summary(g2f_plantdate)
head(g2f_plantdate)

g2f_latlongdf <- merge(g2f_plantdate, g2f_loc1)
g2f_latlongdf$Date.Planted <- as.Date(g2f_latlongdf$Date.Planted, "%m/%d/%y")
g2f_latlongdf$Date.Harvested <- as.Date(g2f_latlongdf$Date.Harvested, "%m/%d/%y")

head(g2f_latlongdf,8)

####NASA power
nasapower_data1 <- data.frame()

for(i in 1:dim(g2f_latlongdf)[1]){
  
  
Year = g2f_latlongdf$Year[i]
LocYear = g2f_latlongdf[i,1]
PD = g2f_latlongdf[i,2]
HD = g2f_latlongdf[i,3]

nasapower_data <- get_power(
    community= "ag",
    pars = c("QV2M",	"T2MDEW",	"PS",	"RH2M",	"WS2M",	"GWETTOP",	"ALLSKY_SFC_SW_DWN",	"ALLSKY_SFC_PAR_TOT",
           "T2M_MAX",	"T2M_MIN",	"T2MWET",	"GWETROOT",	"T2M",	"GWETPROF",	"ALLSKY_SFC_SW_DNI",	"PRECTOTCORR"),
    lonlat = c(g2f_latlongdf[i,7], g2f_latlongdf[i,6]),
    temporal_api = "daily",
    #ates = c(paste0(Year,"-01-01"), paste0(Year,"-12-31")),
    dates = c(PD, HD)
)

nasapower_data <- cbind(LocYear, PD, HD, nasapower_data)
nasapower_data1 <- rbind(nasapower_data1, nasapower_data)

print(i)
}

length(unique(nasapower_data1$LocYear))

summary(nasapower_data1)
nasapower_data1$GDD  <- "NA"

###Calculate GDD
for(i in 1:dim(nasapower_data1)[1]){
  if(nasapower_data1$T2M_MAX[i]>30 & nasapower_data1$T2M_MIN[i]>30){
    nasapower_data1$GDD[i] =  ((30+30)/2)-10 
  } else if(nasapower_data1$T2M_MAX[i]>30 & nasapower_data1$T2M_MIN[i] <10){
    nasapower_data1$GDD[i] =  ((30+10)/2)-10 
  } else if(nasapower_data1$T2M_MAX[i]<10 | nasapower_data1$T2M_MIN[i] <10){
    nasapower_data1$GDD[i] =  ((10+10)/2)-10
  } else if(nasapower_data1$T2M_MAX[i]>30){
    nasapower_data1$GDD[i] =  ((30+nasapower_data1$T2M_MIN[i])/2)-10 
  } else if(nasapower_data1$T2M_MIN[i] <10){
    nasapower_data1$GDD[i] =  ((nasapower_data1$T2M_MAX[i]+10)/2)-10 
  } else {
    nasapower_data1$GDD[i] =  ((nasapower_data1$T2M_MAX[i]+nasapower_data1$T2M_MIN[i])/2)-10
  }
}   
    
head(nasapower_data1)
nasapower_data1$GDD <- as.numeric(nasapower_data1$GDD)    
summary(nasapower_data1)    


nasapower_data1$GDD_Accumulation  <- -999

###Calcuate Accumulated GDD
for(i in 1:dim(nasapower_data1)[1]){
    if(i==1){
      nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i]
  } else if(nasapower_data1$DOY[i] - nasapower_data1$DOY[i-1]==1){
    nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i] + nasapower_data1$GDD_Accumulation[i-1]
  } else {
    nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i]
  }
}
summary(nasapower_data1)

nasapower_data1$GDD_Group = -999

for(i in 1:dim(nasapower_data1)[1]){
  if(nasapower_data1$GDD_Accumulation[i] < 300){
    nasapower_data1$GDD_Group[i] = "<300"
  } else if(nasapower_data1$GDD_Accumulation[i] < 600){
    nasapower_data1$GDD_Group[i] = "300-600"
  } else if(nasapower_data1$GDD_Accumulation[i] < 900){
    nasapower_data1$GDD_Group[i] = "600-900"
  } else if(nasapower_data1$GDD_Accumulation[i] < 1200){
    nasapower_data1$GDD_Group[i] = "900-1200"
  } else if(nasapower_data1$GDD_Accumulation[i] < 1500){
    nasapower_data1$GDD_Group[i] = "1200-1500"
  } else {
    nasapower_data1$GDD_Group[i] = ">1500"
  }
}



for(i in 1:dim(nasapower_data1)[1]){
  if(nasapower_data1$GDD_Accumulation[i] < 500){
    nasapower_data1$GDD_Group[i] = "<500"
  } else if(nasapower_data1$GDD_Accumulation[i] < 1000){
    nasapower_data1$GDD_Group[i] = "500-1000"
  } else {
    nasapower_data1$GDD_Group[i] = ">1000"
  }
}


summary(nasapower_data1)

nasapower_data1$val = -999
for(i in 1:dim(nasapower_data1)[1]){
    if(nasapower_data1$LocYear[i]=="ONH2_2017"){
        nasapower_data1$val[i] = 1
    } else if(nasapower_data1$LocYear[i]=="KSH1_2015"){
      nasapower_data1$val[i] = 2
    } else if(nasapower_data1$LocYear[i]=="MOH1_2015"){
      nasapower_data1$val[i] = 3
    } else if(nasapower_data1$LocYear[i]=="GAH2_2016"){
      nasapower_data1$val[i] = 4
    } else {
      nasapower_data1$val[i] = 5
    }
}
               
ns <- subset(nasapower_data1, nasapower_data1$LocYear=="ONH2_2017")
max(ns$GDD_Accumulation)

col=c("royalblue4", "firebrick", "royalblue4", "firebrick", "coral", "royalblue4", "royalblue4", "green")


library(scales)
n2 <- nasapower_data1 %>%
  arrange(HD)
head(n2)
ggplot(data=nasapower_data1, aes(x=DOY, y=GDD_Accumulation, col=as.factor(val), alpha=as.factor(val))) +
  geom_point() + 
  theme_bw() +
  scale_color_manual(values=c("red", "coral", "green", "blue", "black"),
                     labels=c("ONH2_2017", "KSH1_2015", "MOH1_2015", "GAH2_2016", "Other")) +
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(legend.box.background = element_rect(colour = "black")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, .2)) +
  guides(alpha="none") +
  labs(x="Day in Year",
       y="GDD Accumulation (C)",
       col="LocYears")
  
Trait="T2M_MAX"
  ggplot(data=nasapower_data1, aes(x=GDD_Accumulation, y=T2M_MAX, col=as.factor(val), alpha=as.factor(val))) +
    geom_point() + 
    theme_bw() +
    scale_color_manual(values=c("red", "coral", "green", "blue", "black"),
                       labels=c("ONH2_2017", "KSH1_2015", "MOH1_2015", "GAH2_2016", "Other")) +
    scale_alpha_manual(values = c(1, 1, 1, 1, .2)) +
    guides(alpha="none") +
    theme(axis.title = element_text(size=16, face="bold")) +
    theme(axis.text.x=element_text(size=16)) +
    theme(axis.text.y=element_text(size=16)) +
    theme(legend.box.background = element_rect(colour = "black")) +
    labs(x="GDD Accumulation",
         y="Max Daily Temperature (C)",
         col="LocYears")
    

  
  
  
  
####Should be done from April to October instead of Planting/Harvesting Since we wouldn't know the planting date
  
  ####NASA power
  nasapower_data1 <- data.frame()
  
  for(i in 1:dim(g2f_latlongdf)[1]){
    
    
    Year = g2f_latlongdf$Year[i]
    LocYear = g2f_latlongdf[i,1]
    PD = g2f_latlongdf[i,2]
    HD = g2f_latlongdf[i,3]
    
    nasapower_data <- get_power(
      community= "ag",
      pars = c("QV2M",	"T2MDEW",	"PS",	"RH2M",	"WS2M",	"GWETTOP",	"ALLSKY_SFC_SW_DWN",	"ALLSKY_SFC_PAR_TOT",
               "T2M_MAX",	"T2M_MIN",	"T2MWET",	"GWETROOT",	"T2M",	"GWETPROF",	"ALLSKY_SFC_SW_DNI",	"PRECTOTCORR"),
      lonlat = c(g2f_latlongdf[i,7], g2f_latlongdf[i,6]),
      temporal_api = "daily",
      #dates = c(PD, HD),
      dates = c(paste0(Year,"-04-01"), paste0(Year,"-10-31"))
    )
    
    nasapower_data <- cbind(LocYear, PD, HD, nasapower_data)
    nasapower_data1 <- rbind(nasapower_data1, nasapower_data)
    
    print(i)
  }
  
  length(unique(nasapower_data1$LocYear))
  
  summary(nasapower_data1)
  
  nasapower_data1$GDD  <- "NA"
  
  ###Calculate GDD
  for(i in 1:dim(nasapower_data1)[1]){
    if(nasapower_data1$T2M_MAX[i]>30 & nasapower_data1$T2M_MIN[i]>30){
      nasapower_data1$GDD[i] =  ((30+30)/2)-10 
    } else if(nasapower_data1$T2M_MAX[i]>30 & nasapower_data1$T2M_MIN[i] <10){
      nasapower_data1$GDD[i] =  ((30+10)/2)-10 
    } else if(nasapower_data1$T2M_MAX[i]<10 | nasapower_data1$T2M_MIN[i] <10){
      nasapower_data1$GDD[i] =  ((10+10)/2)-10
    } else if(nasapower_data1$T2M_MAX[i]>30){
      nasapower_data1$GDD[i] =  ((30+nasapower_data1$T2M_MIN[i])/2)-10 
    } else if(nasapower_data1$T2M_MIN[i] <10){
      nasapower_data1$GDD[i] =  ((nasapower_data1$T2M_MAX[i]+10)/2)-10 
    } else {
      nasapower_data1$GDD[i] =  ((nasapower_data1$T2M_MAX[i]+nasapower_data1$T2M_MIN[i])/2)-10
    }
  }   
  
  head(nasapower_data1)
  nasapower_data1$GDD <- as.numeric(nasapower_data1$GDD)    
  summary(nasapower_data1)    
  
  
  nasapower_data1$GDD_Accumulation  <- -999
  
  ###Calcuate Accumulated GDD
  for(i in 1:dim(nasapower_data1)[1]){
    if(i==1){
      nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i]
    } else if(nasapower_data1$DOY[i] - nasapower_data1$DOY[i-1]==1){
      nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i] + nasapower_data1$GDD_Accumulation[i-1]
    } else {
      nasapower_data1$GDD_Accumulation[i] = nasapower_data1$GDD[i]
    }
  }
  summary(nasapower_data1)
  
  

g2f_weather <- nasapower_data1 %>%
  group_by(LocYear, MM) %>%
  summarise(sum_ALLSKY_SFC_PAR_TOT = sum(ALLSKY_SFC_PAR_TOT, na.rm=T),
            avg_T2MWET = mean(T2MWET, na.rm=T),
            avg_QV2M = mean(QV2M, na.rm=T),
            avg_RH2M = mean(RH2M, na.rm=T),
            max_T2M_MAX = max(T2M_MAX, na.rm=T),
            avg_PS = mean(PS, na.rm=T),
            #avg_GWETPROF = mean(GWETPROF, na.rm=T),
            min_T2MDEW = min(T2MDEW, na.rm=T),
            #avg_GWETTOP = mean(GWETTOP, na.rm=T),
            max_WS2M = max(WS2M, na.rm=T),
            min_T2M_MIN = min(T2M_MIN, na.rm=T),
            avg_T2M = mean(T2M, na.rm=T),
            #avg_GWETROOT = mean(GWETROOT, na.rm=T),
            sum_PRECTOTCORR = sum(PRECTOTCORR, na.rm=T))

g2f_GDDMAX <- nasapower_data1 %>%
  group_by(LocYear) %>%
  summarise(MaxGDD = max(GDD_Accumulation, na.rm=T))


g2f_weather1 <- g2f_weather %>%
  pivot_wider(names_from = MM, values_from = c(sum_ALLSKY_SFC_PAR_TOT, 
                                                  avg_T2MWET,
                                                  avg_QV2M,
                                                  avg_RH2M,
                                                  max_T2M_MAX,
                                                  avg_PS,
                                                  #avg_GWETPROF,
                                                  min_T2MDEW,
                                                  #avg_GWETTOP,
                                                  max_WS2M,
                                                  min_T2M_MIN,
                                                  avg_T2M,
                                                  #avg_GWETROOT,
                                                  sum_PRECTOTCORR))




head(nasapower_data1)
head(g2f_weather1)

###Adding in soil data from NRCS 
latlongdf <- g2f_latlongdf[,c(1,6,7)]
colnames(latlongdf) <- c("id", "latitude", "longitude")

g2f_soil <- fetchSoilGrids(x=latlongdf,
                           loc.names=c("id", "latitude", "longitude"),
                           verbose=F,
                           progress=F)
g2f_soil1 <- as.data.frame(g2f_soil@horizons)

g2f_soil1[1:10,1:10]
g2f_soil2 <- reshape(g2f_soil1, direction="wide", timevar = "label",idvar="id")
g2f_soil3 <- merge(latlongdf, g2f_soil2, by="id")
g2f_soil3[1:10,1:10]
colnames(g2f_soil3)[1] <- "LocYear"

g2f_soil3
zeroVar <- function(data, useNA = 'ifany') {
  out <- apply(data, 2, function(x) {length(table(x, useNA = useNA))})
  which(out==1)
}


g2f_soil4 <- g2f_soil3[,-zeroVar(g2f_soil3)]
colnames(g2f_soil4)
g2f_soil5 <- g2f_soil4[,c(1:3,21,26,36,41,46,
                          66,71,81,86,91,
                          111,116,126,131,
                          136)]
head(g2f_soil5)
head(g2f_weather1)
g2f_weather_soil <- merge(g2f_soil5, g2f_weather1, by="LocYear", all.x=TRUE)


g2f_weather_soil1 <- merge(g2f_GDDMAX, g2f_weather_soil, by="LocYear")


head(pheno)
pheno14 <- subset(pheno, pheno$Year==2014)
pheno15 <- subset(pheno, pheno$Year==2015)
pheno16 <- subset(pheno, pheno$Year==2016)
pheno17 <- subset(pheno, pheno$Year==2017)

length(unique(pheno14$LocYear))
length(unique(pheno15$LocYear))
length(unique(pheno16$LocYear))
length(unique(pheno17$LocYear))

length(unique(g2f_weather_soil$LocYear))
EnvGrad <- pheno[,c(2,3,7)]
EnvGrad <- EnvGrad %>%
  distinct(LocYear, .keep_all = TRUE)

g2f_weather_soil_forPRED <- merge(EnvGrad, g2f_weather_soil1, by="LocYear")


customRF <- list(type = "Regression",
                 library = "randomForest",
                 loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                  class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree)
}

#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes




cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
# train model
train_control <- trainControl(method="repeatedcv", 
                              number=10, 
                              repeats=10,
                              allowParallel = TRUE)


#tunegrid <- expand.grid(.mtry = (1:30)) 
tunegrid <- expand.grid(.mtry=c(1:30),.ntree=c(1, 10, 50, 75, 100, 200, 300, 400, 500, 600, 1000, 1500, 2000, 2500))

PREDICTION_ALL <- data.frame()
cor_Pred <- data.frame()

for(i in 2014:2017){
  training <- subset(g2f_weather_soil_forPRED, g2f_weather_soil_forPRED$Year!=i)
  testing <- subset(g2f_weather_soil_forPRED, g2f_weather_soil_forPRED$Year==i)
  
  rownames(training) <- training$LocYear
  rownames(testing) <- testing$LocYear
  
  training <- training[,-c(1:2)]
  testing <- testing[,-c(1:2)]
  
  
  
 # for(ntree in c(1, 10, 50, 75, 100, 200, 300, 400, 500, 600, 1000, 1500, 2000, 2500)){
    set.seed(24)
    g2f_env_pred_train <- train(EnvGrad ~., 
                                data = training,
                                method = customRF,
                                metric = "RMSE",
                                trControl = train_control,  
                                tuneGrid=tunegrid,
                                preProcess = c("center","scale"))
    
  #View the model
  print(g2f_env_pred_train)
  
  plot(varImp(g2f_env_pred_train,scale=T),30)
  
  PredEnvGrad <- predict(g2f_env_pred_train, newdata=testing)
  COR_PREDICTION <- cor(PredEnvGrad, testing$EnvGrad)
  
  cor_int <- cbind(i, COR_PREDICTION)
  cor_Pred <- rbind(cor_Pred, cor_int)
  
  PREDICTION <- as.data.frame(PredEnvGrad)
  PREDICTION_ALL <- rbind(PREDICTION_ALL, PREDICTION)
  
}

cor_Pred
PREDICTION_ALL$LocYear <- rownames(PREDICTION_ALL)
head(PREDICTION_ALL)
test1 <- merge(PREDICTION_ALL, g2f_weather_soil_forPRED, by="LocYear")
cor(test1$PredEnvGrad, test1$EnvGrad)
head(test1)
plot(test1$EnvGrad, test1$PredEnvGrad)


test14 <- subset(test1, test1$Year==2014)
test15 <- subset(test1, test1$Year==2015)
test16 <- subset(test1, test1$Year==2016)
test17 <- subset(test1, test1$Year==2017)

plot(test14$PredEnvGrad, test14$EnvGrad)
plot(test15$PredEnvGrad, test15$EnvGrad)
plot(test16$PredEnvGrad, test16$EnvGrad)
plot(test17$PredEnvGrad, test17$EnvGrad)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###############################[Creating Dataset for Prediction]###################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
pheno <- read.table("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/G2F_pheno.txt", header=T)
head(pheno)

pheno14 <- pheno
for(i in 1:39305){
  if(pheno14[i,3] == 2014){
    pheno14[i,6] = -999
  }
}

pheno15 <- pheno
for(i in 1:39305){
  if(pheno15[i,3] == 2015){
    pheno15[i,6] = -999
  }
}


pheno16 <- pheno
for(i in 1:39305){
  if(pheno16[i,3] == 2016){
    pheno16[i,6] = -999
  }
}



pheno17 <- pheno
for(i in 1:39305){
  if(pheno17[i,3] == 2017){
    pheno17[i,6] = -999
  }
}

write.table(pheno14, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/G2F_pheno14.txt",
            quote=F, sep = " ", row.names = F)

write.table(pheno15, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/G2F_pheno15.txt",
            quote=F, sep = " ", row.names = F)

write.table(pheno16, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/G2F_pheno16.txt",
            quote=F, sep = " ", row.names = F)

write.table(pheno17, "~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/G2F_pheno17.txt",
            quote=F, sep = " ", row.names = F)

length(unique(pheno14$LocYear))
pheno_pred <- merge(pheno, PREDICTION_ALL, by="LocYear")
pheno_pred <- pheno_pred[,-c(8:11)]
length(unique(pheno_pred$LocYear))

pheno14 <- pheno_pred
head(pheno14)
for(i in 1:39305){
  if(pheno14[i,3] == 2014){
    pheno14[i,6] = -999
    pheno14[i,7] = pheno14[i,8]
  }
}

pheno15 <- pheno_pred
for(i in 1:39305){
  if(pheno15[i,3] == 2015){
    pheno15[i,6] = -999
    pheno15[i,7] = pheno15[i,8]
  }
}


pheno16 <- pheno_pred
for(i in 1:39305){
  if(pheno16[i,3] == 2016){
    pheno16[i,6] = -999
    pheno16[i,7] = pheno16[i,8]
  }
}



pheno17 <- pheno_pred
for(i in 1:39305){
  if(pheno17[i,3] == 2017){
    pheno17[i,6] = -999
    pheno17[i,7] = pheno17[i,8]
  }
}



datasets = list(pheno14, pheno15, pheno16, pheno17)
length(datasets)
for(i in 1:4){
  pheno_data <- datasets[[i]]
  
  data.min <- min(pheno_data$EnvGrad,na.rm=TRUE)
  data.max <- max(pheno_data$EnvGrad,na.rm=TRUE)
  print(data.min)
  print(data.max)
  pheno_data$EnvGradStandardized <- -1+2*((pheno_data$EnvGrad-data.min)/(data.max-data.min))
  
  leg4coef <- legendre.polynomials(n=2, normalized=T)
  leg2 <- as.matrix(as.data.frame(polynomial.values(polynomials=leg4coef, x=pheno_data$EnvGradStandardized)))
  
  colnames(leg2) <- c("leg0", "leg1", "leg2")
  leg <- leg2[, 1:ncol(leg2)]
  
  pheno_data <- cbind(pheno_data, leg)
  
  
  
  print(head(pheno_data))
  
  assign(paste0("pheno", 13+i, "_EnvPred"), pheno_data)
  
  write.table(pheno_data, paste0("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/pheno", 13+i, "_EnvPred.txt"),
              quote=F, sep = " ", row.names = F)
}





###Remove genos in testing from training







pheno_pred <- merge(pheno, PREDICTION_ALL, by="LocYear")
pheno_pred <- pheno_pred[,-c(8:11)]
length(unique(pheno_pred$LocYear))

pheno14 <- pheno_pred
head(pheno14)
for(i in 1:39305){
  if(pheno14[i,3] == 2014){
    pheno14[i,6] = -999
    pheno14[i,7] = pheno14[i,8]
  }
}
    
L_Ped <- unique(pheno_pred$PedShort)
pheno14.1 <- data.frame()
for(j in 1:length(L_Ped)){
      PS = L_Ped[j]
      d14 <- subset(pheno14, pheno14$PedShort==PS)

      if(any(d14[,6] < 0)){
        for(k in 1:length(d14[,6])){
          d14[k,6] = -999
        }
      }
    pheno14.1 <- rbind(pheno14.1, d14)  
}
  

pheno14.1 <- pheno14.1 %>%
  arrange(LocYear, PedShort)



pheno15 <- pheno_pred
for(i in 1:39305){
  if(pheno15[i,3] == 2015){
    pheno15[i,6] = -999
    pheno15[i,7] = pheno15[i,8]
  }
}


L_Ped <- unique(pheno_pred$PedShort)
pheno15.1 <- data.frame()
for(j in 1:length(L_Ped)){
  PS = L_Ped[j]
  d15 <- subset(pheno15, pheno15$PedShort==PS)
  
  if(any(d15[,6] < 0)){
    for(k in 1:length(d15[,6])){
      d15[k,6] = -999
    }
  }
  pheno15.1 <- rbind(pheno15.1, d15)  
}


pheno15.1 <- pheno15.1 %>%
  arrange(LocYear, PedShort)



pheno16 <- pheno_pred
for(i in 1:39305){
  if(pheno16[i,3] == 2016){
    pheno16[i,6] = -999
    pheno16[i,7] = pheno16[i,8]
  }
}


L_Ped <- unique(pheno_pred$PedShort)
pheno16.1 <- data.frame()
for(j in 1:length(L_Ped)){
  PS = L_Ped[j]
  d16 <- subset(pheno16, pheno16$PedShort==PS)
  
  if(any(d16[,6] < 0)){
    for(k in 1:length(d16[,6])){
      d16[k,6] = -999
    }
  }
  pheno16.1 <- rbind(pheno16.1, d16)  
}


pheno16.1 <- pheno16.1 %>%
  arrange(LocYear, PedShort)




pheno17 <- pheno_pred
for(i in 1:39305){
  if(pheno17[i,3] == 2017){
    pheno17[i,6] = -999
    pheno17[i,7] = pheno17[i,8]
  }
}


L_Ped <- unique(pheno_pred$PedShort)
pheno17.1 <- data.frame()
for(j in 1:length(L_Ped)){
  PS = L_Ped[j]
  d17 <- subset(pheno17, pheno17$PedShort==PS)
  
  if(any(d17[,6] < 0)){
    for(k in 1:length(d17[,6])){
      d17[k,6] = -999
    }
  }
  pheno17.1 <- rbind(pheno17.1, d17)  
}


pheno17.1 <- pheno17.1 %>%
  arrange(LocYear, PedShort)


datasets = list(pheno14.1, pheno15.1, pheno16.1, pheno17.1)
length(datasets)
for(i in 1:4){
  pheno_data <- datasets[[i]]
  
  data.min <- min(pheno_data$EnvGrad,na.rm=TRUE)
  data.max <- max(pheno_data$EnvGrad,na.rm=TRUE)
  print(data.min)
  print(data.max)
  pheno_data$EnvGradStandardized <- -1+2*((pheno_data$EnvGrad-data.min)/(data.max-data.min))
  
  leg4coef <- legendre.polynomials(n=2, normalized=T)
  leg2 <- as.matrix(as.data.frame(polynomial.values(polynomials=leg4coef, x=pheno_data$EnvGradStandardized)))
  
  colnames(leg2) <- c("leg0", "leg1", "leg2")
  leg <- leg2[, 1:ncol(leg2)]
  
  pheno_data <- cbind(pheno_data, leg)
  
  
  
  print(head(pheno_data))
  
  assign(paste0("pheno", 13+i, "_EnvPred"), pheno_data)
  
  write.table(pheno_data, paste0("~/Documents/Purdue/PhD/Research/Chp4_RRM_EnvironmentalGradient/DATA/forBLUPF90/GenPredPheno/pheno", 13+i, "_EnvPred_NoGenoTraining.txt"),
              quote=F, sep = " ", row.names = F)
}

