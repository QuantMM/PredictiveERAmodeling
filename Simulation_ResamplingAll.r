## -----------------
## Simulations, 2020 March
## -----------------
## Requires R version 3.4. for
library(mlr)
library(MASS)
library(ggplot2)
library(reshape2)

dr <- "C:/Users/SUNMEE/..."
setwd(dr)
source("functions.r")

## -----------------
## Paramters considered: np and startN
## -----------------

Niter <- 500
seedid <- sample.int(100000,Niter)

np <- c(2, 4, 6, 8)						# of predictors per component: for model complexity
startN <- c(50,500)				# sample sizes
senarios <- expand.grid(Nsample = startN, NumX = np)

ResultsCVs <- vector("list",nrow(senarios))
RESULT_LOO <- RESULT_BOOT50 <- RESULT_BOOT20 <- vector("list", nrow(senarios))
RESULT_632_20 <- RESULT_632 <- vector("list", nrow(senarios))

start.time <- Sys.time()
#for (s in 1:nrow(senarios)){
for (s in 1:nrow(senarios)){

	npred <- senarios[s,"NumX"]
	NN <- senarios[s,"Nsample"]

	# Space for saving results
	RESULT_CV3 <- RESULT_CV5 <- RESULT_CV10 <- matrix(, nrow=Niter, ncol=4)
		colnames(RESULT_CV3) <- colnames(RESULT_CV5) <- colnames(RESULT_CV10) <- c("f0","f1","f2","f3")
	ResultsCVs[[s]] <- list(RESULT_CV3, RESULT_CV5, RESULT_CV10)
	
	RESULT_LOO[[s]] <- matrix(, nrow=Niter, ncol=4)
		colnames(RESULT_LOO[[s]]) <- c("f0","f1","f2","f3")
	
	RESULT_BOOT50[[s]] <- RESULT_BOOT20[[s]] <- matrix(, nrow=Niter, ncol=4)
		colnames(RESULT_BOOT20[[s]]) <- colnames(RESULT_BOOT50[[s]]) <- c("f0","f1","f2","f3")
	
	RESULT_632[[s]] <- RESULT_632_20[[s]] <- matrix(, nrow=Niter, ncol=4)
		colnames(RESULT_632[[s]]) <- colnames(RESULT_632_20[[s]]) <- c("f0","f1","f2","f3")
	
	for(i in 1:Niter){
	
		# control seed for sampling in each iteration
		set.seed(seedid[i])
		
		## ------
		## Part1. ERA_tr Data generation by: np && startN
		## ------
		
		# Generated dataset based on ERA structure:
		COVgenerate <- ERAgen(npred = npred)
		Entire <- MASS::mvrnorm(n = NN, mu = rep(0, nrow(COVgenerate$COV)), Sigma = COVgenerate$COV)
			Entire <- data.frame(Entire)
			colnames(Entire)[ncol(Entire)] <- "Y"	#last column = Y
		
		# Predictor sets for F1, F2, ..., F5
		F1Xs <- Entire[,c(paste0("X", 1:npred)),drop=FALSE]
		F2Xs <- Entire[,c(paste0("X", (npred+1):(npred+npred))),drop=FALSE]
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
		## Part3. Resampling strategies: holdout, CV (3, 5, 10), RepCV, LOO, Bootstrap
		## ------

		## fit1 -> True fit: Y = F1 + F2
		X_fit1 <- cbind(F1Xs,F2Xs)
		entireY <- Entire[,"Y",drop=F]
		nvar_fit1 <- matrix( c(ncol(F1Xs),ncol(F2Xs)), ncol=2 )
		
			# scaled data for RMSE
			SCALED_fit1 <- ERA_ts(y=Entire[,"Y",drop=F], X=X_fit1, nvar=nvar_fit1)
			scaledY_fit1 <- SCALED_fit1$y
			scaledX_fit1 <- SCALED_fit1$X
		
		## fit0 -> Underfit: Y = F1 only
		X_fit0 <- F1Xs
		nvar_fit0 <- matrix( c(ncol(F1Xs)), ncol=1 )
		
			# scaled data for RMSE
			SCALED_fit0 <- ERA_ts(y=Entire[,"Y",drop=F], X=X_fit0, nvar=nvar_fit0)
			scaledY_fit0 <- SCALED_fit0$y
			scaledX_fit0 <- SCALED_fit0$X
		
		## fit2 -> Overfit1: Y = F1 + F2 + F1F2
		X_fit2 <- cbind(F1Xs,F2Xs,F3Xs)
		nvar_fit2 <- matrix( c(ncol(F1Xs),ncol(F2Xs),ncol(F3Xs)), ncol=3 )
		
			# scaled data for RMSE
			SCALED_fit2 <- ERA_ts(y=Entire[,"Y",drop=F], X=X_fit2, nvar=nvar_fit0)
			scaledY_fit2 <- SCALED_fit2$y
			scaledX_fit2 <- SCALED_fit2$X
		
		## fit3 -> Overfit2: Y = F1 + F2 + F1F2 + F1^2 + F2^2
		X_fit3 <- cbind(F1Xs,F2Xs,F3Xs,F4Xs,F5Xs)
		nvar_fit3 <- matrix( c(ncol(F1Xs),ncol(F2Xs),ncol(F3Xs),ncol(F4Xs),ncol(F5Xs)), ncol=5 )
		
			# scaled data for RMSE
			SCALED_fit3 <- ERA_ts(y=Entire[,"Y",drop=F], X=X_fit3, nvar=nvar_fit0)
			scaledY_fit3 <- SCALED_fit3$y
			scaledX_fit3 <- SCALED_fit3$X
			
		# ----------
		# N-fold-CV
		# ----------
		CVoptions <- c(3,5,10)
		for (NF in 1:3) {
		
			nfold <- CVoptions[NF]
			CV <- makeResampleDesc("CV", iters = nfold)
			CVid <- makeResampleInstance(CV, size = nrow(Entire))
	
			Err1 <- Err0 <- Err2 <- Err3 <- c()
			for (nf in 1:nfold){
			
				trid <- CVid$train.inds[[nf]]
				tsid <- CVid$test.inds[[nf]]
				
				fit1_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit1[trid,], nvar=nvar_fit1))
				fit0_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit0[trid,], nvar=nvar_fit0))
				fit2_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit2[trid,], nvar=nvar_fit2))
				fit3_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit3[trid,], nvar=nvar_fit3))
				
				if(!is(fit1_tr, "try-error")){				
				ErrTs1 <- RMSE(scaledX_fit1[tsid,], scaledY_fit1[tsid], fit1_tr$W, fit1_tr$A)}
				
				if(!is(fit0_tr, "try-error")){
				ErrTs0 <- RMSE(scaledX_fit0[tsid,], scaledY_fit0[tsid], fit0_tr$W, fit0_tr$A)}
				
				if(!is(fit2_tr, "try-error")){
				ErrTs2 <- RMSE(scaledX_fit2[tsid,], scaledY_fit2[tsid], fit2_tr$W, fit2_tr$A)}
				
				if(!is(fit3_tr, "try-error")){
				ErrTs3 <- RMSE(scaledX_fit3[tsid,], scaledY_fit3[tsid], fit3_tr$W, fit3_tr$A)}
				
				Err1 <- c(Err1,ErrTs1)
				Err0 <- c(Err0,ErrTs0)
				Err2 <- c(Err2,ErrTs2)
				Err3 <- c(Err3,ErrTs3)
				
				rm(fit1_tr);rm(fit0_tr);rm(fit2_tr);rm(fit3_tr)	
			}
			ResultsCVs[[s]][[NF]][i,"f1"] <- mean(Err1)
			ResultsCVs[[s]][[NF]][i,"f0"] <- mean(Err0)
			ResultsCVs[[s]][[NF]][i,"f2"] <- mean(Err2)
			ResultsCVs[[s]][[NF]][i,"f3"] <- mean(Err3)
		}
		rm(trid);rm(tsid);rm(Err1);rm(Err0);rm(Err2);rm(Err3);rm(nf)
	
		# ----------
		# LOO
		# ----------
		LOO <- makeResampleDesc("LOO")
		LOOid <- makeResampleInstance(LOO, size = nrow(Entire))
		
		Err1 <- Err0 <- Err2 <- Err3 <- c()
		for (nf in 1:NN){
			
			trid <- LOOid$train.inds[[nf]]
			tsid <- LOOid$test.inds[[nf]]
			
			fit1_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit1[trid,], nvar=nvar_fit1))
			fit0_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit0[trid,], nvar=nvar_fit0))
			fit2_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit2[trid,], nvar=nvar_fit2))
			fit3_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit3[trid,], nvar=nvar_fit3))
				
			if(!is(fit1_tr, "try-error")){				
			ErrTs1 <- RMSE(scaledX_fit1[tsid,,drop=F], scaledY_fit1[tsid], fit1_tr$W, fit1_tr$A)}
				
			if(!is(fit0_tr, "try-error")){
			ErrTs0 <- RMSE(scaledX_fit0[tsid,,drop=F], scaledY_fit0[tsid], fit0_tr$W, fit0_tr$A)}
				
			if(!is(fit2_tr, "try-error")){
			ErrTs2 <- RMSE(scaledX_fit2[tsid,,drop=F], scaledY_fit2[tsid], fit2_tr$W, fit2_tr$A)}
				
			if(!is(fit3_tr, "try-error")){
			ErrTs3 <- RMSE(scaledX_fit3[tsid,,drop=F], scaledY_fit3[tsid], fit3_tr$W, fit3_tr$A)}
				
			Err1 <- c(Err1,ErrTs1)
			Err0 <- c(Err0,ErrTs0)
			Err2 <- c(Err2,ErrTs2)
			Err3 <- c(Err3,ErrTs3)
				
			rm(fit1_tr);rm(fit0_tr);rm(fit2_tr);rm(fit3_tr)	
		}
		RESULT_LOO[[s]][i,"f1"] <- mean(Err1)
		RESULT_LOO[[s]][i,"f0"] <- mean(Err0)
		RESULT_LOO[[s]][i,"f2"] <- mean(Err2)
		RESULT_LOO[[s]][i,"f3"] <- mean(Err3)
		
		rm(trid);rm(tsid);rm(Err1);rm(Err0);rm(Err2);rm(Err3);rm(nf)
		
		# ----------
		# Bootstrap (iteration20) and .632 estimator
		# ----------
		BOO20 <- makeResampleDesc("Bootstrap", iters = 20)
		BOO20id <- makeResampleInstance(BOO20, size = nrow(Entire))

		Err1 <- Err0 <- Err2 <- Err3 <- c()
		Err1_632 <- Err0_632 <- Err2_632 <- Err3_632 <- c()
		for (nf in 1:20){
		
			trid <- BOO20id$train.inds[[nf]]
			tsid <- BOO20id$test.inds[[nf]]
			
			fit1_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit1[trid,], nvar=nvar_fit1))
			fit0_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit0[trid,], nvar=nvar_fit0))
			fit2_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit2[trid,], nvar=nvar_fit2))
			fit3_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit3[trid,], nvar=nvar_fit3))
				
			if(!is(fit1_tr, "try-error")){				
			ErrTs1 <- RMSE(scaledX_fit1[tsid,,drop=F], scaledY_fit1[tsid], fit1_tr$W, fit1_tr$A)
			ErrTs1_632 <- (.368)*(RMSE(scaledX_fit1[trid,,drop=F], scaledY_fit1[trid], fit1_tr$W, fit1_tr$A))+((.632)*ErrTs1)
			}
				
			if(!is(fit0_tr, "try-error")){
			ErrTs0 <- RMSE(scaledX_fit0[tsid,,drop=F], scaledY_fit0[tsid], fit0_tr$W, fit0_tr$A)
			ErrTs0_632 <- (.368)*(RMSE(scaledX_fit0[trid,,drop=F], scaledY_fit0[trid], fit0_tr$W, fit0_tr$A))+((.632)*ErrTs0)
			}
				
			if(!is(fit2_tr, "try-error")){
			ErrTs2 <- RMSE(scaledX_fit2[tsid,,drop=F], scaledY_fit2[tsid], fit2_tr$W, fit2_tr$A)
			ErrTs2_632 <- (.368)*(RMSE(scaledX_fit2[trid,,drop=F], scaledY_fit2[trid], fit2_tr$W, fit2_tr$A))+((.632)*ErrTs2)
			}
				
			if(!is(fit3_tr, "try-error")){
			ErrTs3 <- RMSE(scaledX_fit3[tsid,,drop=F], scaledY_fit3[tsid], fit3_tr$W, fit3_tr$A)
			ErrTs3_632 <- (.368)*(RMSE(scaledX_fit3[trid,,drop=F], scaledY_fit3[trid], fit3_tr$W, fit3_tr$A))+((.632)*ErrTs3)
			}
				
			Err1 <- c(Err1,ErrTs1)
			Err0 <- c(Err0,ErrTs0)
			Err2 <- c(Err2,ErrTs2)
			Err3 <- c(Err3,ErrTs3)
			
			Err1_632 <- c(Err1_632,ErrTs1_632)
			Err0_632 <- c(Err0_632,ErrTs0_632)
			Err2_632 <- c(Err2_632,ErrTs2_632)
			Err3_632 <- c(Err3_632,ErrTs3_632)
				
			rm(fit1_tr);rm(fit0_tr);rm(fit2_tr);rm(fit3_tr)	
		}
		RESULT_BOOT20[[s]][i,"f1"] <- mean(Err1)
		RESULT_BOOT20[[s]][i,"f0"] <- mean(Err0)
		RESULT_BOOT20[[s]][i,"f2"] <- mean(Err2)
		RESULT_BOOT20[[s]][i,"f3"] <- mean(Err3)
		
		RESULT_632_20[[s]][i,"f1"] <- mean(Err1_632)
		RESULT_632_20[[s]][i,"f0"] <- mean(Err0_632)
		RESULT_632_20[[s]][i,"f2"] <- mean(Err2_632)
		RESULT_632_20[[s]][i,"f3"] <- mean(Err3_632)
		
		rm(trid);rm(tsid);rm(Err1);rm(Err0);rm(Err2);rm(Err3);rm(nf);rm(Err1_632); rm(Err0_632); rm(Err2_632); rm(Err3_632)
		
		# ----------
		# Bootstrap (iteration50) and .632 estimator
		# ----------
		BOO50 <- makeResampleDesc("Bootstrap", iters = 50)
		BOO50id <- makeResampleInstance(BOO50, size = nrow(Entire))
		
		Err1 <- Err0 <- Err2 <- Err3 <- c()
		Err1_632 <- Err0_632 <- Err2_632 <- Err3_632 <- c()
		for (nf in 1:50){
		
			trid <- BOO50id$train.inds[[nf]]
			tsid <- BOO50id$test.inds[[nf]]
			
			fit1_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit1[trid,], nvar=nvar_fit1))
			fit0_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit0[trid,], nvar=nvar_fit0))
			fit2_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit2[trid,], nvar=nvar_fit2))
			fit3_tr <- try(ERA_tr(y=entireY[trid,,drop=F], X=X_fit3[trid,], nvar=nvar_fit3))
				
			if(!is(fit1_tr, "try-error")){				
			ErrTs1 <- RMSE(scaledX_fit1[tsid,,drop=F], scaledY_fit1[tsid], fit1_tr$W, fit1_tr$A)
			ErrTs1_632 <- (.368)*(RMSE(scaledX_fit1[trid,,drop=F], scaledY_fit1[trid], fit1_tr$W, fit1_tr$A))+((.632)*ErrTs1)
			}
				
			if(!is(fit0_tr, "try-error")){
			ErrTs0 <- RMSE(scaledX_fit0[tsid,,drop=F], scaledY_fit0[tsid], fit0_tr$W, fit0_tr$A)
			ErrTs0_632 <- (.368)*(RMSE(scaledX_fit0[trid,,drop=F], scaledY_fit0[trid], fit0_tr$W, fit0_tr$A))+((.632)*ErrTs0)
			}
				
			if(!is(fit2_tr, "try-error")){
			ErrTs2 <- RMSE(scaledX_fit2[tsid,,drop=F], scaledY_fit2[tsid], fit2_tr$W, fit2_tr$A)
			ErrTs2_632 <- (.368)*(RMSE(scaledX_fit2[trid,,drop=F], scaledY_fit2[trid], fit2_tr$W, fit2_tr$A))+((.632)*ErrTs2)
			}
				
			if(!is(fit3_tr, "try-error")){
			ErrTs3 <- RMSE(scaledX_fit3[tsid,,drop=F], scaledY_fit3[tsid], fit3_tr$W, fit3_tr$A)
			ErrTs3_632 <- (.368)*(RMSE(scaledX_fit3[trid,,drop=F], scaledY_fit3[trid], fit3_tr$W, fit3_tr$A))+((.632)*ErrTs3)
			}
				
			Err1 <- c(Err1,ErrTs1)
			Err0 <- c(Err0,ErrTs0)
			Err2 <- c(Err2,ErrTs2)
			Err3 <- c(Err3,ErrTs3)
			
			Err1_632 <- c(Err1_632,ErrTs1_632)
			Err0_632 <- c(Err0_632,ErrTs0_632)
			Err2_632 <- c(Err2_632,ErrTs2_632)
			Err3_632 <- c(Err3_632,ErrTs3_632)
				
			rm(fit1_tr);rm(fit0_tr);rm(fit2_tr);rm(fit3_tr)	
		}
		RESULT_BOOT50[[s]][i,"f1"] <- mean(Err1)
		RESULT_BOOT50[[s]][i,"f0"] <- mean(Err0)
		RESULT_BOOT50[[s]][i,"f2"] <- mean(Err2)
		RESULT_BOOT50[[s]][i,"f3"] <- mean(Err3)
		
		RESULT_632[[s]][i,"f1"] <- mean(Err1_632)
		RESULT_632[[s]][i,"f0"] <- mean(Err0_632)
		RESULT_632[[s]][i,"f2"] <- mean(Err2_632)
		RESULT_632[[s]][i,"f3"] <- mean(Err3_632)
		
		rm(trid);rm(tsid);rm(Err1);rm(Err0);rm(Err2);rm(Err3);rm(nf);rm(Err1_632); rm(Err0_632); rm(Err2_632); rm(Err3_632)
	
	save.image("Resampling1000_from1to5_new.RData")
	} # each iteration ends here

}
end.time <- Sys.time()
(time.taken <- end.time - start.time)

