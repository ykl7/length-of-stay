library(philipsampclient)
library(rjson)
library(ggplot2)
library(caret)
library(MASS)
library(party)
library(stats)
library(rpart)
#initialization api for AMP-ML library
InitializeAMPClientConfig()
outputid_list <- AMPTableMatch(time="NULL", inDatabaseID1="patientID", inDatabaseID2='patientID',
                               locationID1='locationID', locationID2='locationID',
                               dataid='Synthetic_adminClaimsDB_FinalDistribution.csv,Synthetic_ICUclinicalDB_FinalDistribution.csv',
                               PrimaryFtrnames='age,gender,race,year,mortality,LOS',
                               errTolerance='1,0,0,0,0,1', enoughMatch=4, outputdataid="out")

# READ DATA

mi.data <- read.csv("Data.csv")

# CHECK VARIABLE NAMES

names(mi.data)

# REMOVE ROWS WITH N/A LOS

los.available <- c(!is.na(mi.data$Los))
available.data <- mi.data[los.available, ]

available.data[complete.cases(available.data), ]

# FIND OPTIMAL NUMBER OF CLUSTERS (ELBOW IN GRAPH)

set.seed(15)
wss <- (nrow(available.data)-1)*sum(apply(available.data, 2, var))
for (i in 2:15) wss[i] <- sum(kmeans(available.data$Los, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# k-MEANS CLUSTERING

los.cluster <- kmeans(available.data$Los, 5, nstart = 500, iter.max = 50000)
los.cluster$cluster <- as.factor(los.cluster$cluster)
ggplot(available.data, aes(available.data$SL_NO, available.data$Los, color = los.cluster$cluster)) + geom_point() + xlab("ID") + ylab("LOS")
available.data$los.group <- los.cluster$centers[los.cluster$cluster]
los.cluster$centers

non.outlier.data <- available.data

non.outlier.data$db.ratio <- non.outlier.data$SBP/non.outlier.data$DBP
max.fbs <- max(non.outlier.data$FBS)
max.bmi <- max(non.outlier.data$BMI)
max.vldl <- max(non.outlier.data$VLDL)
max.frs <- max(non.outlier.data$FRS)
non.outlier.data$FRS <- non.outlier.data$FRS/max.frs
non.outlier.data$BMI <- non.outlier.data$BMI/max.bmi
non.outlier.data$VLDL <- non.outlier.data$VLDL/max.vldl
non.outlier.data$FBS <- non.outlier.data$FBS/max.fbs
non.outlier.data$NT.proBNP <- as.factor(non.outlier.data$NT.proBNP)
non.outlier.data$commonality <- (non.outlier.data$Anxiety.Depression + non.outlier.data$Pain.Discomfort - non.outlier.data$Mobility)

non.outlier.data$Hosp_No <- non.outlier.data$D.O.A <- non.outlier.data$X <- non.outlier.data$X.1 <- non.outlier.data$X.2 <- non.outlier.data$Angina <- 
  non.outlier.data$HTN <- non.outlier.data$DM <- non.outlier.data$Dyslipidemia <- non.outlier.data$vessel.involved <- non.outlier.data$vessel.involved.1 <- 
  non.outlier.data$Diagnosis <- non.outlier.data$Smoking <- non.outlier.data$Tobacco <- non.outlier.data$Alcohol <- non.outlier.data$Gutka <- 
  non.outlier.data$family.h.o <- non.outlier.data$Pre...Physical.Ex <- non.outlier.data$Gly.HB <- non.outlier.data$pre.Bpcheck.up <- 
  non.outlier.data$pre.BS.check.up <- non.outlier.data$Usual.activitis <- non.outlier.data$Self.care <- non.outlier.data$HA <-
  non.outlier.data$ACS_Score <- non.outlier.data$NT.Pro.BNP <- non.outlier.data$Diet <- non.outlier.data$pre.diet <- non.outlier.data$Education <- NULL

# RANDOM FOREST

output.forest.part <- rpart(non.outlier.data$Los ~ non.outlier.data$Age + non.outlier.data$SES + non.outlier.data$profession + non.outlier.data$Income + 
                         non.outlier.data$w.H.ratio + non.outlier.data$HTN.DUR + non.outlier.data$DM.DUR + non.outlier.data$DYSP.DUR, data = non.outlier.data, method = "anova")
plot(output.forest.part)
text(output.forest.part)
printcp(output.forest.part)
plotcp(output.forest.part)
print(output.forest.part)
rsq.rpart(output.forest.part)

rForest <- rpart (non.outlier.data$Los ~ non.outlier.data$Angina.Class + non.outlier.data$TC.HDLratio + non.outlier.data$db.ratio + non.outlier.data$TG + 
                    non.outlier.data$FBS, data = non.outlier.data, method = "anova")
plot(rForest)
text(rForest)
printcp(rForest)
plotcp(rForest)
print(rForest)
rsq.rpart(rForest)

# DATASET SPLIT FUNCTION

splitdataframe <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/2))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

partitioned.data <- splitdataframe(non.outlier.data, seed = 15)
test.data <- partitioned.data$testset
train.data <- partitioned.data$trainset

dim(train.data)

los.predict.lm <- lm(Los ~ commonality + db.ratio + Angina.Class + 
                       TC.HDLratio + (profession * Income * Age) + 
                       (Trop..T * NT.proBNP) + FBS + TG + 
                       BMI + EF + FRS, data = train.data)
summary(los.predict.lm)
stepAIC(los.predict.lm, steps = 1000, direction = "both")
plot(los.predict.lm)
