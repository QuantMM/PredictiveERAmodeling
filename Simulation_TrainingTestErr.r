install.packages("MASS")
library(MASS)
## -----------------
## Simulations, 2020 March
## -----------------

dr <- "C:/Users/SUNMEE/..."
setwd(dr)
source("functions.r")

## -----------------
## Paramters considered: np and startN
## -----------------

Niter <- 5000
seedid <- sample.int(100000,Niter)

np <- c(2, 4, 6, 8)						# of predictors per component: for model complexity
startN <- c(50,100,200,500,1000)		# sample sizes
senarios <- expand.grid(NumX = np, Nsample = startN)

RESULT_TrainingErr <- vector("list", nrow(senarios))
RESULT_TestErr1 <- vector("list", nrow(senarios))		# test set size = Nsample
RESULT_TestErr2 <- vector("list", nrow(senarios))		# test set size = 10 * Nsample

start.time <- Sys.time()
for (s in 1:nrow(senarios)){

	npred <- senarios[s,"NumX"]
	NN <- senarios[s,"Nsample"]

	# Space for saving results
	RESULT_TrainingErr[[s]]$ad.R <- matrix(, nrow=Niter, ncol=4)
	RESULT_TrainingErr[[s]]$RMSE <- matrix(, nrow=Niter, ncol=4) # ave diff b/t observedYs & values predicted by the model
		colnames(RESULT_TrainingErr[[s]]$ad.R) <- colnames(RESULT_TrainingErr[[s]]$RMSE) <- c("f0","f1","f2","f3")
	
	for(i in 1:Niter){
	
		# control seed for sampling in each iteration
		set.seed(seedid[i])
		
		## ------
		## Part1. ERA_tr Data generation by: np && startN
		## ------
		
		# Generated train/test datasets based on ERA structure:
		COVgenerate <- ERAgen(npred = npred)
		
		# Test set set-ups
		TEST1 <- MASS::mvrnorm(n = NN, mu = rep(0, nrow(COVgenerate$COV)), Sigma = COVgenerate$COV)
			TEST1 <- data.frame(TEST1); colnames(TEST1)[ncol(TEST1)] <- "Y"
			F1F2Xs_ts1 <-  subset(TEST1, select = -Y)
		TEST2 <- MASS::mvrnorm(n = NN*10, mu = rep(0, nrow(COVgenerate$COV)), Sigma = COVgenerate$COV)
			TEST2 <- data.frame(TEST2); colnames(TEST2)[ncol(TEST2)] <- "Y"
			F1F2Xs_ts2 <- subset(TEST2, select = -Y)
		
		# Train set set-ups
		TRAIN <- MASS::mvrnorm(n = NN, mu = rep(0, nrow(COVgenerate$COV)), Sigma = COVgenerate$COV)
			TRAIN <- data.frame(TRAIN)
			colnames(TRAIN)[ncol(TRAIN)] <- "Y"	#last column = Y
		
		# Predictor sets for F1, F2, ..., F5
		F1Xs <- TRAIN[,c(paste0("X", 1:npred)),drop=FALSE]
		F2Xs <- TRAIN[,c(paste0("X", (npred+1):(npred+npred))),drop=FALSE]
			# Generate interaction terms b/w Xs in F1 and F2:
			xpairs <- expand.grid(colnames(F1Xs),colnames(F2Xs), stringsAsFactors=F)
			inters <- data.frame(matrix(, nrow= nrow(F1Xs), ncol= nrow(xpairs)))
				for (tt in 1:nrow(xpairs)) {
				
				factor_1 <- xpairs[tt, 1]
				factor_2 <- xpairs[tt, 2]
				inters[,tt] <- F1Xs[,colnames(F1Xs) == factor_1] * F2Xs[,colnames(F2Xs) == factor_2]
				colnames(inters)[tt] <- gsub("\\.", ":", paste(factor_1, factor_2, sep=":"))
				}
		F3Xs <- inters
		F4Xs <- '^'(F1Xs,2)
		F5Xs <- '^'(F2Xs,2)
	
		## ------
		## Part2. Model evaluation -- training error (for qunatifying model performance)
		## : calc adjustedR^2 and AIC as in traditional ERA
		## : calc RMSEs of four models on the same data the model was trained on
		## ------
		
		entireY <- TRAIN[,"Y",drop=F]
		test1Y <- TEST1[,"Y",drop=F]
		test2Y <- TEST2[,"Y",drop=F]
				
		## fit1 -> True fit: Y = F1 + F2
		nFs <- 2
		X_fit1 <- cbind(F1Xs,F2Xs)
		nvar_fit1 <- matrix( c(ncol(F1Xs),ncol(F2Xs)), ncol=nFs )
		fit1 <- try(ERA_tr(y=entireY, X=X_fit1, nvar=nvar_fit1))
		
		if(!is(fit1, "try-error")){
			#results:
			
			# TrainErr
			scaledY <- ERA_ts(y=entireY, X=X_fit1, nvar=nvar_fit1)$y
			scaledX <- fit1$X
			RESULT_TrainingErr[[s]]$RMSE[i,"f1"] <- RMSE(scaledX, scaledY, fit1$W, fit1$A)
			
			# TestErr based on fit1
			scaled_test1Y <- ERA_ts(y=test1Y, X=F1F2Xs_ts1, nvar=nvar_fit1)$y
			scaled_test1X <- ERA_ts(y=test1Y, X=F1F2Xs_ts1, nvar=nvar_fit1)$X
			scaled_test2Y <- ERA_ts(y=test2Y, X=F1F2Xs_ts2, nvar=nvar_fit1)$y
			scaled_test2X <- ERA_ts(y=test2Y, X=F1F2Xs_ts2, nvar=nvar_fit1)$X
			RESULT_TestErr1[[s]] <- c(RESULT_TestErr1[[s]], RMSE(scaled_test1X, scaled_test1Y, fit1$W, fit1$A))
			RESULT_TestErr2[[s]] <- c(RESULT_TestErr2[[s]], RMSE(scaled_test2X, scaled_test2Y, fit1$W, fit1$A))
			
			# R_ad
			fit1Data <- data.frame(Y=fit1$adj.DV, X=fit1$F);reg1 <- lm(Y~.-1,data=fit1Data)
			p1 <- sum(nvar_fit1)+nFs # how many parameters in fit1
			RESULT_TrainingErr[[s]]$ad.R[i,"f1"] <- 1-{ (1-summary(reg1)$r.squared)*((NN-1)*(1/(NN-p1-1)))}
		}
		
		## fit0 -> Underfit: Y = F1 only
		nFs <- 1
		X_fit0 <- F1Xs
		nvar_fit0 <- matrix( c(ncol(F1Xs)), ncol=nFs )
		fit0 <- try(ERA_tr(y=entireY, X=F1Xs, nvar=nvar_fit0))
		
		if(!is(fit0, "try-error")){
			#results:
			scaledX <- fit0$X
			RESULT_TrainingErr[[s]]$RMSE[i,"f0"] <- RMSE(scaledX, scaledY, fit0$W, fit0$A)
			
			fit0Data <- data.frame(Y=fit0$adj.DV, X=fit0$F); reg0 <- lm(Y~.-1,data=fit0Data)
			p0 <- sum(nvar_fit0)+nFs
			RESULT_TrainingErr[[s]]$ad.R[i,"f0"] <- 1-{ (1-summary(reg0)$r.squared)*((NN-1)*(1/(NN-p0-1)))}
		}
		
		## fit2 -> Overfit1: Y = F1 + F2 + F1F2
		nFs <- 3
		X_fit2 <- cbind(F1Xs,F2Xs,F3Xs)
		nvar_fit2 <- matrix( c(ncol(F1Xs),ncol(F2Xs),ncol(F3Xs)), ncol=nFs )
		fit2 <- try(ERA_tr(y=entireY, X=X_fit2, nvar=nvar_fit2))
		
		if(!is(fit2, "try-error")){
			#results:
			scaledX <- fit2$X
			RESULT_TrainingErr[[s]]$RMSE[i,"f2"] <- RMSE(scaledX, scaledY, fit2$W, fit2$A)
			
			fit2Data <- data.frame(Y=fit2$adj.DV, X=fit2$F); reg2 <- lm(Y~.-1,data=fit2Data)
			p2 <- sum(nvar_fit2)+nFs
			RESULT_TrainingErr[[s]]$ad.R[i,"f2"] <- 1-{ (1-summary(reg2)$r.squared)*((NN-1)*(1/(NN-p2-1)))}
		}
		
		## fit3 -> Overfit2: Y = F1 + F2 + F1F2 + F1^2 + F2^2
		nFs <- 5
		X_fit3 <- cbind(F1Xs,F2Xs,F3Xs,F4Xs,F5Xs)
		nvar_fit3 <- matrix( c(ncol(F1Xs),ncol(F2Xs),ncol(F3Xs),ncol(F4Xs),ncol(F5Xs)), ncol=nFs )
		fit3 <- try(ERA_tr(y=entireY, X=X_fit3, nvar=nvar_fit3))

		if(!is(fit3, "try-error")){
			#results:
			scaledX <- fit3$X
			RESULT_TrainingErr[[s]]$RMSE[i,"f3"] <- RMSE(scaledX, scaledY, fit3$W, fit3$A)
			
			fit3Data <- data.frame(Y=fit3$adj.DV, X=fit3$F); reg3 <- lm(Y~.-1,data=fit3Data)
			p3 <- sum(nvar_fit3)+nFs
			RESULT_TrainingErr[[s]]$ad.R[i,"f3"] <- 1-{ (1-summary(reg3)$r.squared)*((NN-1)*(1/(NN-p3-1)))}
		}
	
	# Remove fits in each iteration
	rm(fit0); rm(fit1); rm(fit2); rm(fit3)
	
	save.image("TrainTestErr.RData")
	} # each iteration ends here
}
end.time <- Sys.time()
(time.taken <- end.time - start.time)

