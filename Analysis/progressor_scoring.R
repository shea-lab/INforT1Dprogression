#Per mouse scoring on -7 versus -5 signature

dmall = cbind.data.frame(dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10)
mtIDs = colnames(dmall)

genes_EN = NODt.5t1sig
batchcorr_counts = t.5t1
t.5t1sig.85 = g_EN.85

ENdm1 <- na.omit(dm1[t.5t1sig.85, ]) 
dm1samp = colnames(ENdm1)
dmtimes = c("t0","t0.5","t1","t2","t3","t4")

ENdm2 <- na.omit(dm2[t.5t1sig.85, ]) 
dm2samp = colnames(ENdm2)
ENdm3 <- na.omit(dm3[t.5t1sig.85, ]) 
dm3samp = colnames(ENdm3)
ENdm4 <- na.omit(dm4[t.5t1sig.85, ]) 
dm4samp = colnames(ENdm4)
ENdm5 <- na.omit(dm5[t.5t1sig.85, ]) 
dm5samp = colnames(ENdm5)
ENdm6 <- na.omit(dm6[t.5t1sig.85, ]) 
dm6samp = colnames(ENdm6)
ENdm7 <- na.omit(dm7[t.5t1sig.85, ]) 
dm7samp = colnames(ENdm7)
ENdm8 <- na.omit(dm8[t.5t1sig.85, ]) 
dm8samp = colnames(ENdm8)
ENdm9 <- na.omit(dm9[t.5t1sig.85, ]) 
dm9samp = colnames(ENdm9)
ENdm10 <- na.omit(dm10[t.5t1sig.85, ]) 
dm10samp = colnames(ENdm10)


timereorder = c(1,6,2,3,4,5)
ENdm1 = ENdm1[, timereorder]
ENdm2 = ENdm2[, timereorder]
ENdm3 = ENdm3[, timereorder]
ENdm4 = ENdm4[, timereorder]

ENdmall = cbind.data.frame(ENdm1,ENdm2,ENdm3,ENdm4,ENdm5,ENdm6,ENdm7,ENdm8,ENdm9,ENdm10)
dmallsamp = colnames(ENdmall)
alltimes = c(rep(c("t0","t0.5","t1","t2","t3","t4"),4),"t0","t1","t2","t3","t4",
	rep(c("t0","t0.5","t1","t2","t3","t4"),2),"t0","t1","t2","t3","t4",
	rep(c("t0","t0.5","t1","t2","t3","t4"),2))


#SVD
n_cent <- ENdmall - rowMeans(ENdmall) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
# Z <- t(n_centered) %*% svd1$v # Z=XV
#svd2$u # U are unscaled PCs
par(mfrow=c(1,1)) # Only one plot per page
plot(svd2$u[,1], svd2$u[,2], col = coul, main = "Genes v Samples (SVD)", pch=19,
     xlab = "V1", ylab = "V2")
text(svd2$u[,1], svd2$u[,2], dmallsamp, cex=0.7, pos=4)
centr_healthy <- colMeans(rbind(svd2$u[alltimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:58]
euc <- cbind.data.frame(dmallsamp, euclid_dist)
rownames(euc) <- c(dmallsamp)
euc_SVD <- cbind(euc, alltimes)

#Colors
colSide <- alltimes
for (x in 1:length(alltimes)) {
  if (colSide[x] == "t0") {colSide[x] <- "skyblue"
  } else if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t2") {colSide[x] <- "mediumvioletred"
  } else if (colSide[x] == "t3") {colSide[x] <- "orangered"
  } else if (colSide[x] == "t4") {colSide[x] <- "red4"
  }}


ggplot(euc, aes(x = reorder(alltimes, euclid_dist), y = euclid_dist, color = names(colSide))) + 
  geom_boxplot(lwd=1.5, fatten=0.75) + xlab("") + ylab("Singular Value Decomposition") + 
  scale_color_manual(name="Group", values=c("darkred", "blue", "orange", "black")) + 
  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

#RF
predictor_data <- t(ENdmall) # Transpose data & assign genes as predictors.
target <- alltimes # Set variable to predict target (reject status)
target[target=="t0"] <- "Early"
target[target=="t0.5"] <- "Early"
target[target=="t1"] <- "Late"
target[target=="t2"] <- "Late"
target[target=="t3"] <- "Late"
target[target=="t4"] <- "Late"
target <- as.factor(target)
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=F)[1]]
sampsizes <- rep(min_size,num_classes)
set.seed(123)
rf_output <- randomForest(x=predictor_data, y=target, importance=T, ntree=10001, 
                          proximity=T, sampsize=sampsizes, na.action=na.omit)
rf_importances <- importance(rf_output, scale=F) # Importance for each variable
confusion <- rf_output$confusion # Calculates sensitivity, specificity, accuracy
sensitivity <- (confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity <- (confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error <- rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy <- 1-overall_error
class1_error <- paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error <- paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy <- 100-overall_error
sens_out <- paste("sensitivity=",sensitivity, sep="")
spec_out <- paste("specificity=",specificity, sep="")
err_out <- paste("overall error rate=",overall_error,sep="")
acc_out <- paste("overall accuracy=",overall_accuracy,sep="")
misclass_1 <- paste(confusion[1,2], rownames(confusion)[1],"misclassified as", 
                    colnames(confusion)[2], sep=" ")
misclass_2 <- paste(confusion[2,1], rownames(confusion)[2],"misclassified as", 
                    colnames(confusion)[1], sep=" ")
confusion_out <- confusion[1:2,1:2]
confusion_out <- cbind(rownames(confusion_out), confusion_out)
p_predictors <- varImpPlot(rf_output, type=2, n.var=12, scale=F, # Top variables
                           main="Variable Importance (Gini) for EN predictors")
target_labels <- as.vector(target) # MDS Class Separation
#target_labels[target_labels=="t0"] <- "0"
target_labels[target_labels=="Early"] <- "E"
target_labels[target_labels=="Late"] <- "L"
#target_labels[target_labels=="Diabetic"] <- "D"
plot_MDS <- MDSplot(rf_output, target,k=2,xlab="",ylab="",pch=target_labels, 
    palette=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), main="MDS plot")
plot_MDS <- plot_MDS$points
plot_MDS <- as.data.frame(plot_MDS)
p_mds <- ggplot(plot_MDS, aes(x=`Dim 1`,y=`Dim 2`, color=target_labels)) + 
  geom_point(aes(shape=donor),size=3) + 
  geom_text(aes(label = time),nudge_x=0.03, nudge_y=-0.01)
centr_healthy <- colMeans(rbind(plot_MDS[target == "Early",]),dims=1)
euclid <- rbind(centr_healthy, plot_MDS) # Euclidean distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:58]
euc <- cbind.data.frame(dmallsamp, euclid_dist)
rownames(euc) <- c(dmallsamp)
euc_RF <- cbind(euc, alltimes)
#ggplot(euc, aes(x = simple, y = euclid_dist, fill = cohort)) + 
#  geom_boxplot() + xlab("") + ylab("") + ggtitle("Random Forest")
ggplot(euc, aes(x = reorder(alltimes, euclid_dist), y = euclid_dist, color = names(colSide))) + 
  geom_boxplot(lwd=1.5, fatten=0.75) + xlab("") + ylab("Random Forest") + theme_classic() + 
  scale_color_manual(name="Group", values=c("darkred", "blue", "orange", "black")) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))
predictions <- as.vector(rf_output$votes[,2]) # ROC Curve, D:H votes as prediction

#ROCR currently supports only evaluation of binary classification tasks.
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
options(digits=2) # Vote Distributions
#out <- histbackback(split(rf_output$votes[,"Diabetic"], target), probability=F, 
#                   axes=T, xlim=c(-50,50), main='Vote distributions for mice classified by RF', 
#                   ylab="Fraction votes (Diabetic)")
#barplot(-out$left, col="red" , horiz=T, space=0, add=T, axes=F)
#barplot(out$right, col="blue", horiz=T, space=0, add=T, axes=F)

# 13 (New) Scoring System ----
par(mar = c(4,4,4,4), oma = c(1,1,1,1))
plot(euc_SVD$euclid_dist, euc_RF$euclid_dist)
scores <- cbind.data.frame(euc_SVD$euclid_dist, euc_RF$euclid_dist)
colnames(scores) <- c("SVD", "RF")
rownames(scores) <- dmallsamp
write.csv(scores, "allsamplesscores.csv", row.names=TRUE)

#Need to add clin_info struct to get this graphing to run
ggplot(scores, aes(x=SVD,y=RF, color=colSide)) + 
  geom_text(aes(label = alltimes),nudge_x=0.03, nudge_y=-0.01) + 
  xlab("Singular Value Decomposition (arb. units)") + 
  ylab("Random Forest prediction, Diabetic probability") + 
  coord_equal()

scores$SVD <- (scores$SVD - min(scores$SVD))/(max(scores$SVD) - min(scores$SVD))
scores$RF <- (scores$RF - min(scores$RF))/(max(scores$RF) - min(scores$RF))

NODRF = scores$RF #1-5 healthy, 6-10 diseased

# All mice and times
par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores$SVD, scores$RF, col=colSide, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.83, .45, legend=c("6 wks old","-7 wks","-5 wks","-3 wks","-1 wk","Onset"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)

# Individual mice

colIn <- dmtimes
for (x in 1:length(dmtimes)) {
  if (colIn[x] == "t0") {colIn[x] <- "skyblue"
  } else if (colIn[x] == "t0.5") {colIn[x] <- "dodgerblue"
  } else if (colIn[x] == "t1") {colIn[x] <- "darkorchid"
  } else if (colIn[x] == "t2") {colIn[x] <- "mediumvioletred"
  } else if (colIn[x] == "t3") {colIn[x] <- "orangered"
  } else if (colIn[x] == "t4") {colIn[x] <- "red4"
  }}


m1 = scores[1:6,]
m2 = scores[7:12,]
m3 = scores[13:18,]
m4 = scores[19:24,]
m5 = scores[25:29,]
m6 = scores[30:35,]
m7 = scores[36:41,]
m8 = scores[42:46,] 
m9 = scores[47:52,]
m10 = scores[53:58,]

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(m9$SVD, m9$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 9")

#Individual mouse SVD scores
n_cent <- ENdm1 - rowMeans(ENdm1) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm1samp, euclid_dist)
rownames(euc) <- c(dm1samp)
euc_SVD1 <- cbind(euc, dmtimes)
scores1 <- cbind.data.frame(euc_SVD1$euclid_dist, euc_RF[1:6,]$euclid_dist)
colnames(scores1) <- c("SVD", "RF")
scores1$SVD <- (scores1$SVD - min(scores1$SVD))/(max(scores1$SVD) - min(scores1$SVD))
scores1$RF <- (scores1$RF - min(scores1$RF))/(max(scores1$RF) - min(scores1$RF))
rownames(scores1)=dm1samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores1$SVD, scores1$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 1")

n_cent <- ENdm2 - rowMeans(ENdm2) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm2samp, euclid_dist)
rownames(euc) <- c(dm2samp)
euc_SVD2 <- cbind(euc, dmtimes)
scores2 <- cbind.data.frame(euc_SVD2$euclid_dist, euc_RF[7:12,]$euclid_dist)
colnames(scores2) <- c("SVD", "RF")
scores2$SVD <- (scores2$SVD - min(scores2$SVD))/(max(scores2$SVD) - min(scores2$SVD))
scores2$RF <- (scores2$RF - min(scores2$RF))/(max(scores2$RF) - min(scores2$RF))
rownames(scores2)=dm2samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores2$SVD, scores2$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 2")

n_cent <- ENdm3 - rowMeans(ENdm3) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm3samp, euclid_dist)
rownames(euc) <- c(dm3samp)
euc_SVD3 <- cbind(euc, dmtimes)
scores3 <- cbind.data.frame(euc_SVD3$euclid_dist, euc_RF[13:18,]$euclid_dist)
colnames(scores3) <- c("SVD", "RF")
scores3$SVD <- (scores3$SVD - min(scores3$SVD))/(max(scores3$SVD) - min(scores3$SVD))
scores3$RF <- (scores3$RF - min(scores3$RF))/(max(scores3$RF) - min(scores3$RF))
rownames(scores3)=dm3samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores3$SVD, scores3$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 3")

n_cent <- ENdm4 - rowMeans(ENdm4) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm4samp, euclid_dist)
rownames(euc) <- c(dm4samp)
euc_SVD4 <- cbind(euc, dmtimes)
scores4 <- cbind.data.frame(euc_SVD4$euclid_dist, euc_RF[19:24,]$euclid_dist)
colnames(scores4) <- c("SVD", "RF")
scores4$SVD <- (scores4$SVD - min(scores4$SVD))/(max(scores4$SVD) - min(scores4$SVD))
scores4$RF <- (scores4$RF - min(scores4$RF))/(max(scores4$RF) - min(scores4$RF))
rownames(scores4)=dm4samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores4$SVD, scores4$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 4")

n_cent <- ENdm6 - rowMeans(ENdm6) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm6samp, euclid_dist)
rownames(euc) <- c(dm6samp)
euc_SVD6 <- cbind(euc, dmtimes)
scores6 <- cbind.data.frame(euc_SVD6$euclid_dist, euc_RF[30:35,]$euclid_dist)
colnames(scores6) <- c("SVD", "RF")
scores6$SVD <- (scores6$SVD - min(scores6$SVD))/(max(scores6$SVD) - min(scores6$SVD))
scores6$RF <- (scores6$RF - min(scores6$RF))/(max(scores6$RF) - min(scores6$RF))
rownames(scores6)=dm6samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores6$SVD, scores6$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 6")

n_cent <- ENdm7 - rowMeans(ENdm7) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm7samp, euclid_dist)
rownames(euc) <- c(dm7samp)
euc_SVD7 <- cbind(euc, dmtimes)
scores7 <- cbind.data.frame(euc_SVD7$euclid_dist, euc_RF[36:41,]$euclid_dist)
colnames(scores7) <- c("SVD", "RF")
scores7$SVD <- (scores7$SVD - min(scores7$SVD))/(max(scores7$SVD) - min(scores7$SVD))
scores7$RF <- (scores7$RF - min(scores7$RF))/(max(scores7$RF) - min(scores7$RF))
rownames(scores7)=dm7samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores7$SVD, scores7$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 7")

n_cent <- ENdm9 - rowMeans(ENdm9) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm9samp, euclid_dist)
rownames(euc) <- c(dm9samp)
euc_SVD9 <- cbind(euc, dmtimes)
scores9 <- cbind.data.frame(euc_SVD9$euclid_dist, euc_RF[47:52,]$euclid_dist)
colnames(scores9) <- c("SVD", "RF")
scores9$SVD <- (scores9$SVD - min(scores9$SVD))/(max(scores9$SVD) - min(scores9$SVD))
scores9$RF <- (scores9$RF - min(scores9$RF))/(max(scores9$RF) - min(scores9$RF))
rownames(scores9)=dm9samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores9$SVD, scores9$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 9")

n_cent <- ENdm10 - rowMeans(ENdm10) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[dmtimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:6]
euc <- cbind.data.frame(dm10samp, euclid_dist)
rownames(euc) <- c(dm10samp)
euc_SVD10 <- cbind(euc, dmtimes)
scores10 <- cbind.data.frame(euc_SVD10$euclid_dist, euc_RF[53:58,]$euclid_dist)
colnames(scores10) <- c("SVD", "RF")
scores10$SVD <- (scores10$SVD - min(scores10$SVD))/(max(scores10$SVD) - min(scores10$SVD))
scores10$RF <- (scores10$RF - min(scores10$RF))/(max(scores10$RF) - min(scores10$RF))
rownames(scores10)=dm10samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores10$SVD, scores10$RF, col=colIn, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t0.5","t1","t2","t3","t4"),
       fill=c("skyblue","dodgerblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 10")


#No t0.5
alttimes = c("t0","t1","t2","t3","t4")
colIn2 <- alttimes
for (x in 1:length(alttimes)) {
  if (colIn2[x] == "t0") {colIn2[x] <- "skyblue"
  } else if (colIn2[x] == "t1") {colIn2[x] <- "darkorchid"
  } else if (colIn2[x] == "t2") {colIn2[x] <- "mediumvioletred"
  } else if (colIn2[x] == "t3") {colIn2[x] <- "orangered"
  } else if (colIn2[x] == "t4") {colIn2[x] <- "red4"
  }}

n_cent <- ENdm5 - rowMeans(ENdm5) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[alttimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:5]
euc <- cbind.data.frame(dm5samp, euclid_dist)
rownames(euc) <- c(dm5samp)
euc_SVD5 <- cbind(euc, alttimes)
scores5 <- cbind.data.frame(euc_SVD5$euclid_dist, euc_RF[25:29,]$euclid_dist)
colnames(scores5) <- c("SVD", "RF")
scores5$SVD <- (scores5$SVD - min(scores5$SVD))/(max(scores5$SVD) - min(scores5$SVD))
scores5$RF <- (scores5$RF - min(scores5$RF))/(max(scores5$RF) - min(scores5$RF))
rownames(scores5)=dm5samp

par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores5$SVD, scores5$RF, col=colIn2, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t1","t2","t3","t4"),
       fill=c("skyblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 5")

n_cent <- ENdm8 - rowMeans(ENdm8) # First, center data on genes
svd2 <- svd(t(n_cent)) # Apply SVD on transposed data
centr_healthy <- colMeans(rbind(svd2$u[alttimes == "t0",]),dims=1)
euclid <- rbind(centr_healthy, svd2$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:5]
euc <- cbind.data.frame(dm8samp, euclid_dist)
rownames(euc) <- c(dm8samp)
euc_SVD8 <- cbind(euc, alttimes)
scores8 <- cbind.data.frame(euc_SVD8$euclid_dist, euc_RF[42:46,]$euclid_dist)
colnames(scores8) <- c("SVD", "RF")
scores8$SVD <- (scores8$SVD - min(scores8$SVD))/(max(scores8$SVD) - min(scores8$SVD))
scores8$RF <- (scores8$RF - min(scores8$RF))/(max(scores8$RF) - min(scores8$RF))
rownames(scores8)=dm8samp

SVDind1 = cbind.data.frame(scores1$SVD,scores2$SVD,scores3$SVD,scores4$SVD,
	scores6$SVD,scores7$SVD,scores9$SVD,scores10$SVD)
rownames(SVDind1) = dmtimes
write.csv(SVDind1, "SVDscoresno58.csv", row.names=TRUE)

SVDind2 = cbind.data.frame(scores5$SVD,scores8$SVD)
rownames(SVDind2) = alttimes
write.csv(SVDind2, "SVDscores58.csv", row.names=TRUE)


par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0.5,0.5,0.5,0.5)) # Adjust margins
plot(scores8$SVD, scores8$RF, col=colIn2, pch=19, cex=2, 
     mgp=c(2.5,1,0), cex.lab=1.25,
     xaxs="i", yaxs="i",
     xlim=c(0,1.2), ylim=c(0,1.1),
     ylab="Random Forest prediction (prob.)", 
     xlab="Singular Value Decomposition (arb. units)") 
legend(.9, .45, legend=c("t0","t1","t2","t3","t4"),
       fill=c("skyblue","darkorchid","mediumvioletred","orangered","red4"), cex=1.3)
title("Mouse 8")


#Training/test of t0.5 vs all later timepoints
deltas = cbind.data.frame(dD.5,dD1,dD2,dD3,dD4)
deltimes= c(rep("t0.5: -7 weeks",length(dD.5[1,])),rep("t1: -5 weeks",length(dD1[1,])),
	rep("t2: -3 weeks",length(dD2[1,])),rep("t3: -1 week",length(dD3[1,])),
	rep("t4: Diabetes onset",length(dD4[1,])))
timegroups = c(rep("-7 weeks",length(dD.5[1,])),rep("late",length(dD1[1,])),
	rep("late",length(dD2[1,])),rep("late",length(dD3[1,])),rep("late",length(dD4[1,])))
grouped7vall = factor(timegroups, labels = c(0,1)) 
allv7samps = deltimes
metadata = cbind(allv7samps, grouped7vall)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)

#Figure 6D, ROC of top 3 alphas for 50 iterations
allv7flexi <- flexiDEG.function4(deltas, metadata, validation_option = 2)

EN7v.1 <- na.omit(deltas[unique(rownames(deltas)[as_vector(allv7flexi[[1]])]), ]) 
EN7v.2 <- na.omit(deltas[unique(rownames(deltas)[as_vector(allv7flexi[[2]])]), ]) 
EN7v.3 <- na.omit(deltas[unique(rownames(deltas)[as_vector(allv7flexi[[3]])]), ]) 

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid","mediumvioletred","orangered","red4")
colSide <- deltimes
for (x in 1:length(deltimes)) {
  if (colSide[x] == "t0.5: -7 weeks") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1: -5 weeks") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t2: -3 weeks") {colSide[x] <- "mediumvioletred"
  } else if (colSide[x] == "t3: -1 week") {colSide[x] <- "orangered"
  } else if (colSide[x] == "t4: Diabetes onset") {colSide[x] <- "red4"
  }}
colind = colSide
names(colind) = deltimes

heatmap.2(as.matrix(EN7v.2), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN7v.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)
