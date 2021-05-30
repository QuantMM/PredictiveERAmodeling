round(m4_t_RMSE,6)
round(m4_v_RMSE,6)
round(m4_cv,6)
round(m4_BOO,6)
round(m4_632,6)

# 4. LOOCV err (using training set)

	NN <- sample(1:nrow(X_t), 6000)
	LOO <- makeResampleDesc("LOO")
	LOOid <- makeResampleInstance(LOO, size = nrow(X_t))
	
	Err4 <- c() # err of model4
	for (nf in 1:length(NN)){
		
		nnf <- NN[nf]	
		trid <- LOOid$train.inds[[nnf]]
		tsid <- LOOid$test.inds[[nnf]]
		
		fit4_tr <- try(ERA_tr(y=y_t[trid,,drop=F], X=X_t[trid,], nvar=nvar_t))
		
		if(!is(fit4_tr, "try-error")){				
		ErrTs4 <- RMSE(X_t[tsid,,drop=F], y_t[tsid,,drop=F], fit4_tr$W, fit4_tr$A)}
		Err4 <- c(Err4,ErrTs4)
		rm(fit4_tr)	
	}
	
	m4_LOOCV <- mean(Err4)
round(m4_LOOCV,6)

