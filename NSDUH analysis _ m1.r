## -----------------
## Requires R version 3.4. for
## -----------------
library(mlr)
library(MASS)
#library(ggplot2)
#library(reshape2)

setwd("C:/...")
D <- read.csv("merged2016-2019.csv", header = T, sep = ",", na.strings=c("NA"))

# ---------------------------------------------------------
# Train (learning set) / Test (validation set)
# ---------------------------------------------------------
D_t <- D[D$year!=2019,]
D_v <- D[D$year==2019,]

# try
D_t <- D_t[sample(1:nrow(D_t),1000),]
D_v <- D_v[sample(1:nrow(D_v),100),]

##########-------------------------------------------------
##########
# Model 1: Y = F1 + F2 
##########
##########-------------------------------------------------

# -----------------
# 1. Normalize Data
# -----------------

# 1-1. Train set:

	y_t <- D_t[,"ndssansp",drop=FALSE]
	x1_t <- D_t[,c("ircigage","iralcage","irmjage"),drop=FALSE]
	x2_t <- D_t[,c("WSPDSC2","WHODASC3","MHSUTK_U","AMDEY2_U"),drop=FALSE]
	X_t <- cbind(x1_t,x2_t)

	nvar_t <- matrix(c(ncol(x1_t),ncol(x2_t)), ncol=2)

	NORMD = cbind(y_t,X_t)
	NORMD <- na.omit(NORMD)
	NORMD_n = nrow(y_t)
	dataN <- scale(NORMD)/sqrt(NORMD_n-1)

	y_t <- dataN[,1,drop=FALSE]
	X_t <- dataN[,-1,drop=FALSE]

# 1-2. Test set:

	y_v <- D_v[,"ndssansp",drop=FALSE]
	x1_v <- D_v[,c("ircigage","iralcage","irmjage"),drop=FALSE]
	x2_v <- D_v[,c("WSPDSC2","WHODASC3","MHSUTK_U","AMDEY2_U"),drop=FALSE]
	X_v <- cbind(x1_v,x2_v)

	nvar_v <- matrix(c(ncol(x1_v),ncol(x2_v)), ncol=2)

	NORMD = cbind(y_v,X_v)
	NORMD <- na.omit(NORMD)
	NORMD_n = nrow(y_v)
	dataN <- scale(NORMD)/sqrt(NORMD_n-1)

	y_v <- dataN[,1,drop=FALSE]
	X_v <- dataN[,-1,drop=FALSE]

# -----------------
# 2. Train / Test set RMSE
# -----------------
source("functions.r")

# 2-1. Train set RMSE
	m4_t <- ERA_tr(y=y_t, X=X_t, nvar=nvar_t)
	m4_t_RMSE <- RMSE(X_t, y_t, m4_t$W, m4_t$A)

# 2-2. Test set RMSE
	m4_v_RMSE <- RMSE(X_v, y_v, m4_t$W, m4_t$A)
		
round(m4_t_RMSE,6);round(m4_v_RMSE,6)


# -----------------
# 3. 10-fold err (using training set)
# -----------------

	nfold <- 10
	CV <- makeResampleDesc("CV", iters = nfold)
	CVid <- makeResampleInstance(CV, size = nrow(X_t))
	
	Err4 <- c() # err of model4
	for (nf in 1:nfold){
	
		trid <- CVid$train.inds[[nf]]
		tsid <- CVid$test.inds[[nf]]
		
		fit4_tr <- try(ERA_tr(y=y_t[trid,,drop=F], X=X_t[trid,], nvar=nvar_t))

		if(!is(fit4_tr, "try-error")){				
		ErrTs4 <- RMSE(X_t[tsid,], y_t[tsid,,drop=F], fit4_tr$W, fit4_tr$A)}
		
		Err4 <- c(Err4,ErrTs4)
		
		rm(fit4_tr)	
	}
	
	m4_cv <- mean(Err4)
	rm(Err4)

# -----------------
# 4. LOOCV err (using training set)
# -----------------

	NN <- nrow(X_t)
	LOO <- makeResampleDesc("LOO")
	LOOid <- makeResampleInstance(LOO, size = nrow(X_t))
	
	Err4 <- c() # err of model4
	for (nf in 1:NN){
		
		trid <- LOOid$train.inds[[nf]]
		tsid <- LOOid$test.inds[[nf]]
		
		fit4_tr <- try(ERA_tr(y=y_t[trid,,drop=F], X=X_t[trid,], nvar=nvar_t))
		
		if(!is(fit4_tr, "try-error")){				
		ErrTs4 <- RMSE(X_t[tsid,,drop=F], y_t[tsid,,drop=F], fit4_tr$W, fit4_tr$A)}
			
		Err4 <- c(Err4,ErrTs4)
		
		rm(fit4_tr)	
	}
	
	m4_LOOCV <- mean(Err4)
	rm(Err4)

# -----------------
# 5. Bootstrap (iteration50): OOB and .632 err (using training set)
# -----------------

	BOO50 <- makeResampleDesc("Bootstrap", iters = 50)
	BOO50id <- makeResampleInstance(BOO50, size = nrow(X_t))
	
	Err4 <- c() # err of model4
	Err4_632 <- c()
	for (nf in 1:50){
	
		trid <- BOO50id$train.inds[[nf]]
		tsid <- BOO50id$test.inds[[nf]]
		
		fit4_tr <- try(ERA_tr(y=y_t[trid,,drop=F], X=X_t[trid,], nvar=nvar_t))

		if(!is(fit4_tr, "try-error")){				
		ErrTs4 <- RMSE(X_t[tsid,,drop=F], y_t[tsid], fit4_tr$W, fit4_tr$A)
		ErrTs4_632 <- (.368)*(RMSE(X_t[trid,,drop=F], y_t[trid], fit4_tr$W, fit4_tr$A))+((.632)*ErrTs4)
		}
			
		Err4 <- c(Err4,ErrTs4)
		Err4_632 <- c(Err4_632,ErrTs4_632)
			
		rm(fit4_tr)
	}
	
	m4_BOO <- mean(Err4)
	rm(Err4)
	
	m4_632 <- mean(Err4_632)
	rm(Err4_632)
	


# ---------------------------------------------------------
# Print outcomes
# ---------------------------------------------------------
round(m4_t_RMSE,6)
round(m4_v_RMSE,6)
m4_cv
m4_LOOCV
m4_BOO
m4_632

##########-------------------------------------------------
##########
# Model 1 FINISHED here
##########
##########-------------------------------------------------
