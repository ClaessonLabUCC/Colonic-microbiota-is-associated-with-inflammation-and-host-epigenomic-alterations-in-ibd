
#####Classification
group <- binary_label

group[,1] <- gsub("3", "1", group[,1])
group[,1] <- gsub("2", "0", group[,1])
label <- as.numeric(group)
x <- Data

#ensuring there are no nas in the Data
#x <- data.matrix(t(na.omit(t(x))))

#adding column names for features so that feature importances can be collected for feature selection
if(length(colnames(x)) == 0 ){colnames(x) <- c(1:ncol(x))}

#making vectors to store data collected from each iteration of the for loop
predictions <- numeric()
importance <- data.frame(as.matrix(colnames(x), nrow = ncol(x), ncol = nrow(x)))
colnames(importance) <- c("names")
predictions <- numeric()
max_depth_list <- numeric()
subsample_list <- numeric()
colsample_bytree_list <- numeric()
min_child_weight_list <- numeric()


########################################################################################################################
#######PART 1, feature selection, Data already includes random noise feature, it is the last feature column#############
########################################################################################################################


#for loop which iterates through each row/sample in the data and leaves it out of training & parameter selection and uses it for testing
for (a in (1:nrow(x))){
  
  
  
  uni <- as.numeric(rownames(subset(Patient_ID, Patient_ID[,1] != Patient_ID[a,1])))
  uni_sel <- numeric()
  for (r in (unique(Patient_ID[uni,1]))) {
    chosen <- sample(rownames(subset(Patient_ID, Patient_ID[,1]== r)), size =1, replace = F)
    uni_sel <- append(uni_sel, chosen, after = length(uni_sel))
  }
  uni_sel <- as.numeric(uni_sel)
  
  
  ######parameter selection loop 
  best_param = list()
  best_seednumber = 1234
  best_auc = -Inf
  best_auc_index = 0
  
  for (iter in 1:10000) {
    seed.number = iter*100
    param <- list(objective = "binary:logistic",
                  eval_metric = "auc",
                  max_depth = sample(10:100, 1),
            subsample = runif(1, .5, .9),
                  colsample_bytree = runif(1, .5, .9)
    )
    
    #the number of cross validation and rounds
    cv.nround = 30
    cv.nfold = 4
    
    set.seed(seed.number)
    mdcv <- xgb.cv(data=x[uni_sel,], params = param, nthread=8, label = label[uni_sel], stratified = TRUE,
                   nfold=cv.nfold, nrounds=cv.nround,
                   verbose = T, early.stopping.round=10, maximize=TRUE)
    
    max_auc = max(mdcv$evaluation_log[, test_auc_mean])
    max_auc_index = which.max(mdcv$evaluation_log[, test_auc_mean])
    
    if (max_auc > best_auc) {
      best_auc = max_auc
      best_auc_index = max_auc_index
      best_seednumber = seed.number
      best_param = param
    }
  }
  
  nround = max_auc_index
  set.seed(best_seednumber)
  
  xgb1 <- xgboost(x[uni_sel,] ,as.numeric(label[uni_sel]),params=best_param, nrounds=nround, nthread=6, verbose = F)
  test <- xgb.dump(xgb1)
  
  write.table(colnames(x), paste0("models/", "fmap",a,".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
  invisible(xgb.dump(model = xgb1, fname = paste0("models/","xgb.dump",a), fmap = paste0("models/","fmap",a,".txt"), with_stats = TRUE))
  
  
  if(test[[2]] == "0:leaf=0"){} else if (test[[2]] == "0:leaf=-0"){} else{tryCatch(
    (imp_dis <- xgb.importance(feature_names = colnames(x), model= xgb1)), silent = FALSE)}              
  imp <- as.matrix(imp_dis$Gain)
  rownames(imp) <- imp_dis$Feature 
  importance <- merge(importance, imp, by.x = "names", by.y = "row.names", all = T)
  test_sample <- xgb.DMatrix( t(data.matrix(x[c(a),])), label = label[[a]])
  predictionsi <- predict(xgb1, test_sample)
  # predictionsi  <- predictionsi <- sum((predictionsi > 0.5) != label[[a]])/length(label[[a]])
  predictions <- append(predictions, predictionsi , after = length(predictions))
  max_depth_list <- append(max_depth_list, param[3] , after = length(max_depth_list))
  subsample_list <- append(subsample_list, param[4] , after = length(subsample_list))
  colsample_bytree_list <- append(colsample_bytree_list, param[5] , after = length(colsample_bytree_list))
  
}




assign(paste("predictions_xgb"),predictions)
pred_label <- cbind(label, predictions_xgb)
#write.csv(pred_label, "predictions.csv")
assign(paste("important_features_xgb"),importance)

###############################################################
####graphing auc from from models without feature selection####
###############################################################

roc <- roc(label, predictions_xgb)
xgb.pred <- prediction(predictions_xgb, label)
xgb.perf <- performance(xgb.pred, "tpr", "fpr")
plot(xgb.perf,
     avg="threshold",
     colorize=TRUE,
     lwd=1,
     print.cutoffs.at=seq(0, 1, by=0.05),
     text.adj=c(-0.5, 0.5),
     text.cex=0.5)
grid(col="lightgray")
axis(1, at=seq(0, 1, by=0.1))
axis(2, at=seq(0, 1, by=0.1))
abline(v=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
abline(h=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
lines(x=c(0, 1), y=c(0, 1), col="black", lty="dotted")

plot(roc(label, predictions_xgb), print.auc=TRUE)
roc2 <- simple_roc(as.numeric(label), predictions_xgb)
with(roc2, lines(1 - FPR, TPR, col="blue", lty=2))
roc
power.roc.test(roc)
roc.area(label, predictions_xgb)



summary(unlist(max_depth_list))
summary(unlist(subsample_list))
summary(unlist(colsample_bytree_list))
summary(unlist(min_child_weight_list))



l <- important_features_xgb
safe <- as.matrix(l[,"names"])
noise_num <- which(l[,"names"] == "noise")
l[is.na(l)] <- 0
for(a in 1:(nrow(l)-1)){for (b in 2:(ncol(l))){if(as.numeric(l[a,b]) > as.numeric(l[noise_num,b])){l[a,b] <- 1} else {l[a,b] <- 0}}}
l <- l[,!(colnames(l)) %in% c("names")] 
l <- as.matrix(apply(l, 2,  as.numeric))
sum_imp <- as.matrix(rowSums(l, na.rm = TRUE))
rownames(sum_imp) <- (1:nrow(sum_imp))
threshold <- sum_imp[noise_num,]
if((is.na(threshold)) == TRUE){threshold = 0}
sum_imp <- subset(sum_imp, sum_imp > threshold)
important <- as.numeric(rownames(sum_imp))
assign(paste("selected_features"),safe[important])
assign(paste("num_appearances"),sum_imp)
store <- cbind( as.matrix(selected_features), num_appearances)

########################################################################################################################
#######PART 2#############
########################################################################################################################
x <- x[,store[,1]]
predictions <- numeric()
importance <- data.frame(as.matrix(colnames(x), nrow = ncol(x), ncol = 1))
colnames(importance) <- c("names")
predictions <- numeric()
max_depth_list <- numeric()
subsample_list <- numeric()
colsample_bytree_list <- numeric()
min_child_weight_list <- numeric()



#for loop which iterates through each row/sample in the data and leaves it out of training & parameter selection and uses it for testing
for (a in (1:nrow(x))){
  
  uni <- as.numeric(rownames(subset(Patient_ID, Patient_ID[,1] != Patient_ID[a,1])))
  uni_sel <- numeric()
  for (r in (unique(Patient_ID[uni,1]))) {
    chosen <- sample(rownames(subset(Patient_ID, Patient_ID[,1]== r)), size =1, replace = F)
    uni_sel <- append(uni_sel, chosen, after = length(uni_sel))
  }
  uni_sel <- as.numeric(uni_sel)
  
  
  ######parameter selection loop 
  best_param = list()
  best_seednumber = 1234
  best_auc = -Inf
  best_auc_index = 0
  
  for (iter in 1:10000) {
    seed.number = iter*100
    param <- list(objective = "binary:logistic",
                  eval_metric = "auc",
                  max_depth = sample(1:100, 1),
                  subsample = runif(1, .5, .8),
                  colsample_bytree = runif(1, .5, .8)
    )
    
    #the number of cross validation and rounds
    cv.nround = 30
    cv.nfold = 2
    
    set.seed(seed.number)
    mdcv <- xgb.cv(data=x[uni_sel,], params = param, nthread=10, label = label[uni_sel], 
                   nfold=cv.nfold, nrounds=cv.nround,
                   verbose = T, early.stopping.round=8, maximize=TRUE, stratified = TRUE)
    
    max_auc = max(mdcv$evaluation_log[, test_auc_mean])
    max_auc_index = which.max(mdcv$evaluation_log[, test_auc_mean])
    
    if (max_auc > best_auc) {
      best_auc = max_auc
      best_auc_index = max_auc_index
      best_seednumber = seed.number
      best_param = param
    }
  }
  
  nround = max_auc_index
  set.seed(best_seednumber)
  
  xgb1 <- xgboost(x[uni_sel,] ,as.numeric(label[uni_sel]),params=best_param, nrounds=nround, nthread=6, verbose = F)
  test <- xgb.dump(xgb1)
  
  write.table(colnames(x), paste0("models/", "fmap",a,".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
  invisible(xgb.dump(model = xgb1, fname = paste0("models/","xgb.dump",a), fmap = paste0("models/","fmap",a,".txt"), with_stats = TRUE))
  
  
  if(test[[2]] == "0:leaf=0"){} else if (test[[2]] == "0:leaf=-0"){} else{tryCatch(
    (imp_dis <- xgb.importance(feature_names = colnames(x), model= xgb1)), silent = FALSE)}              
  imp <- as.matrix(imp_dis$Gain)
  rownames(imp) <- imp_dis$Feature 
  importance <- merge(importance, imp, by.x = "names", by.y = "row.names", all = T)
  test_sample <- xgb.DMatrix( t(data.matrix(x[c(a),])), label = label[[a]])
  predictionsi <- predict(xgb1, test_sample)
  # predictionsi  <- predictionsi <- sum((predictionsi > 0.5) != label[[a]])/length(label[[a]])
  predictions <- append(predictions, predictionsi , after = length(predictions))
  max_depth_list <- append(max_depth_list, param[3] , after = length(max_depth_list))
  subsample_list <- append(subsample_list, param[4] , after = length(subsample_list))
  colsample_bytree_list <- append(colsample_bytree_list, param[5] , after = length(colsample_bytree_list))
  
  
}





assign(paste("predictions_xgb"),predictions)
pred_label <- cbind(label, predictions_xgb)
assign(paste("important_features_xgb"),importance)

###############################################################
####graphing auc from from models without feature selection####
###############################################################

roc <- roc(label, predictions_xgb)
xgb.pred <- prediction(predictions_xgb, label)
xgb.perf <- performance(xgb.pred, "tpr", "fpr")
plot(xgb.perf,
     avg="threshold",
     colorize=TRUE,
     lwd=1,
     print.cutoffs.at=seq(0, 1, by=0.05),
     text.adj=c(-0.5, 0.5),
     text.cex=0.5)
grid(col="lightgray")
axis(1, at=seq(0, 1, by=0.1))
axis(2, at=seq(0, 1, by=0.1))
abline(v=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
abline(h=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
lines(x=c(0, 1), y=c(0, 1), col="black", lty="dotted")

plot(roc(label, predictions_xgb), print.auc=TRUE)
roc2 <- simple_roc(as.numeric(label), predictions_xgb)
with(roc2, lines(1 - FPR, TPR, col="blue", lty=2))
roc
power.roc.test(roc)
roc.area(label, predictions_xgb)



summary(unlist(max_depth_list))
summary(unlist(subsample_list))
summary(unlist(colsample_bytree_list))
summary(unlist(min_child_weight_list))



l <- important_features_xgb
safe <- as.matrix(l[,"names"])
noise_num <- which(l[,"names"] == "noise")
l[is.na(l)] <- 0
for(a in 1:(nrow(l)-1)){for (b in 2:(ncol(l))){if(as.numeric(l[a,b]) > as.numeric(l[noise_num,b])){l[a,b] <- 1} else {l[a,b] <- 0}}}
l <- l[,!(colnames(l)) %in% c("names")] 
l <- as.matrix(apply(l, 2,  as.numeric))
sum_imp <- as.matrix(rowSums(l, na.rm = TRUE))
rownames(sum_imp) <- (1:nrow(sum_imp))
threshold <- sum_imp[noise_num,]
if((is.na(threshold)) == TRUE){threshold = 0}
sum_imp <- subset(sum_imp, sum_imp > threshold)
important <- as.numeric(rownames(sum_imp))
assign(paste("selected_features"),safe[important])
assign(paste("num_appearances"),sum_imp)
store <- cbind( as.matrix(selected_features), num_appearances)
