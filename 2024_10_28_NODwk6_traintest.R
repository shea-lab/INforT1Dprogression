#NOD v NS at all times
#Unnormalized
subNOD = clin_info_NOD$strain_subsets
groups = factor(subNOD, labels = c(0,1)) #0: NS, 1: NOD

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold", "skyblue")
colSide <- subNOD
for (x in 1:length(subNOD)) {
  if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"  
  } else if (colSide[x] == "NOD progressors") {colSide[x] <- "skyblue"
  }}
colind = colSide
names(colind) = subNOD

  xfactors <- cbind.data.frame(groups) 
samplesbc = colnames(NOD_filt)
  rownames(xfactors) <- samplesbc
  dataset <- t(NOD_filt) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 1, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODvNSallt = genes_EN
batchcorr_counts = NOD_filt
pdf(file= "NOD vs NS all times, bc2.pdf" )

#Normalized to t0
NOD_filt_norm = cbind.data.frame(deltas, nsdeltas)
delsubNOD = c(rep("NOD progressors", 48), rep("NOD non-progressors", 22))
groups = factor(delsubNOD, labels = c(0,1)) #0: NS, 1: NOD

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold", "skyblue")
colSide <- delsubNOD
for (x in 1:length(delsubNOD)) {
  if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"  
  } else if (colSide[x] == "NOD progressors") {colSide[x] <- "skyblue"
  }}
colind = colSide
names(colind) = delsubNOD

  xfactors <- cbind.data.frame(groups) 
samplesbc = colnames(NOD_filt_norm)
  rownames(xfactors) <- samplesbc
  dataset <- t(NOD_filt_norm) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 1, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODvNSalltnorm = genes_EN
batchcorr_counts = NOD_filt_norm
pdf(file= "NOD vs NS all times norm, bc2.pdf" )

# Wk 6 only Random Forest for ROC ---- 
genes_EN = NODvNSwk6bc2
g_EN.9wk6 <- unique(rownames(wk6)[genes_EN$ENgenes$`0.9`])
  n_EN.9 <- na.omit(wk6[g_EN.9, ])
predictor_data <- t(n_EN.9) # Transpose data & assign genes as predictors.
target <- clin_info_wk6$strain_subsets # Set variable to predict target (reject status)
target[target=="NOD progressors"] <- "Progressor"
target[target=="NOD non-progressors"] <- "Non-progressor"
target <- as.factor(target)
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=F)[1]]
sampsizes <- rep(min_size,num_classes)
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
target_labels[target_labels=="NOD progressors"] <- "P"
target_labels[target_labels=="NOD non-progressors"] <- "N"
plot_MDS <- MDSplot(rf_output, target,k=2,xlab="",ylab="",pch=target_labels, 
                    palette=c("red", "blue"), main="MDS plot")
plot_MDS <- plot_MDS$points
plot_MDS <- as.data.frame(plot_MDS)
p_mds <- ggplot(plot_MDS, aes(x=`Dim 1`,y=`Dim 2`, color=target_labels)) + 
  geom_point(aes(shape=donor),size=3) + 
  geom_text(aes(label = time),nudge_x=0.03, nudge_y=-0.01)
centr_healthy <- colMeans(rbind(plot_MDS[target == "NOD non-progressors",]),dims=1)
euclid <- rbind(centr_healthy, plot_MDS) # Euclidean distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:21]
euc <- cbind.data.frame(NODwk6samples, euclid_dist)
rownames(euc) <- c(rownames(clin_info_wk6))
euc_RF <- cbind(euc, NODwk6samples)
#ggplot(euc, aes(x = simple, y = euclid_dist, fill = cohort)) + 
#  geom_boxplot() + xlab("") + ylab("") + ggtitle("Random Forest")
ggplot(euc, aes(x = reorder(simple, euclid_dist), y = euclid_dist, color = names(colSide))) + 
  geom_boxplot(lwd=1.5, fatten=0.75) + xlab("") + ylab("Random Forest") + theme_classic() + 
  scale_color_manual(name="Group", values=c("darkred", "blue", "orange", "black")) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))
predictions <- as.vector(rf_output$votes[,2]) # ROC Curve, D:H votes as prediction
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
options(digits=2) # Vote Distributions
out <- histbackback(split(rf_output$votes[,"Diabetic"], target), probability=F, 
                    axes=T, xlim=c(-50,50), main='Vote distributions for mice classified by RF', 
                    ylab="Fraction votes (Diabetic)")
barplot(-out$left, col="red" , horiz=T, space=0, add=T, axes=F)
barplot(out$right, col="blue", horiz=T, space=0, add=T, axes=F)


groups = factor(clin_info_wk6$strain_subsets, labels = c(0,1)) #1: NOD, 0:NS
 xfactors <- cbind.data.frame(groups) 
samplesbc = colnames(wk6)
  rownames(xfactors) <- samplesbc
  dataset <- t(wk6) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 0.6, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable

 cvfit <- suppressMessages(cv.glmnet(x, y, alpha = 0.9, 
                                        nfolds = 100, family = "multinomial"))
    best.lambda <- cvfit$lambda.min
    fit <- suppressMessages(glmnet(x, y, alpha = 0.9, lambda = best.lambda, 
                  family = "multinomial"))
 probabilities <-(predict(fit, x, type = "response", s = best.lambda))
      # Create a binary matrix for the true labels of the test data
      labels <- as.factor(test.data$groups)
      # Flatten the probabilities matrix to match the dimensions of the labels
      prob_flat <- probabilities
   # Convert the probability matrix to a data frame for binomial classification
        prob_flat <- as.data.frame(prob_flat) #Binomial
        # Calculate the probability of the second class for binary classification
        prob_flat$s2<- 1 - prob_flat$s1
        # Rename columns to represent class labels
        colnames(prob_flat) <- c("0", "1")
        # Convert the probability data frame back to a matrix
        prob_flat<-as.matrix(prob_flat)
          
pred <- prediction(prob_flat,groups)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))


#Simple EN, NOD vs NS week 6 all run in signature
groups = factor(clin_info_wk6$strain_subsets, labels = c(0,1)) #1: NOD, 0:NS
  xfactors <- cbind.data.frame(groups) 
samplesbc = colnames(wk6)
  rownames(xfactors) <- samplesbc
  dataset <- t(wk6) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 0.6, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODvNSwk6bc2.6train = genes_EN
 wk6.test = rownames(test.data)
 test = wk6[ , wk6.test]
 batchcorr_counts = test

wk6_sub = clin_info_wk6[wk6.test, ]
colSide <- wk6_sub$strain_subsets
for (x in 1:length(colSide)) {
  if (colSide[x] == "NOD progressors") {colSide[x] <- "skyblue"  
  } else if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"
  }}
colind = colSide
names(colind) = wk6_sub$strain_subsets

pdf(file= "NOD week 6, simple train and test, bc2.pdf" )

