## -------------------
## Functions for simulation
## -------------------

ERAgen <- function(corX=0, b1=0.3, b2=0.5568, npred=npred, r=0){

	b_vec <- matrix(c(b1,b2),nrow=2)
	w1 <- sample(seq(0.8,1,0.01), npred, replace=T)
	w2 <- sample(seq(0.8,1,0.01), npred, replace=T)
	
	# cov of predictors
	cx <- matrix(1,nrow=npred,ncol=npred)*corX - diag(1,npred)*corX + diag(1,npred)
	
	# cross-corr among predictors and Y
	W1 <- cbind(w1,0)
	W2 <- cbind(0,w2)
	W <- rbind(W1,W2)
	k <- 2 # ncomp
	CX <- kronecker(diag(k), cx)
	SdCX <- sqrt(t(W)%*%CX%*%W)
	Ws <- W *(1/SdCX[1])
	W_vec <- Ws[Ws!=0]
	
	crossCor <- CX%*%Ws%*%b_vec

	# pop. cov matrix
	COVX <- cbind(CX, crossCor)
	COVY <- cbind(t(crossCor),1)
	COV <- rbind(COVX,COVY)
	
	out.ERA_gen <- list(b_vec=b_vec,W_vec=W_vec,COV=COV)
	res <- out.ERA_gen
	res$b_vec
	res$W_vec
	res$COV
	invisible(res)
}


ERA_ts <- function(y=y, X=X, nvar=nvar, dist=1){
	
	if(dist==1){ Yfam = "gaussian"
	} else if (dist==2) { Yfam = "binomial" #(link = "logit")
	} else if (dist==3) { Yfam = "poisson" #(link = "log")
	} else if (dist==4) { Yfam = "Gamma" #(link = "inverse")
	} else if (dist==5) { Yfam = "inverse.gaussian" #(link = "1/mu^2")
	}
	
	data = cbind(y,X)
	ndset = length(nvar)
	ncase = nrow(y)
	ny = ncol(y)
	sum_nvar = sum(nvar)
	# ------------------
	# data normalization
	# ------------------
	dataN <- scale(data)/sqrt(ncase-1)
	y <- dataN[,1]
	X <- dataN[,-1]
	
	output <- list(X=X,y=y)
	res <- output
	res$X
	res$y
	
	class(res) <- c("ERA_ts", class(res))
	invisible(res)
}


RMSE <- function(XX=XX, YY=YY, W=W, A=A){

	y_obs <- YY
	y_hat <- as.matrix(XX) %*% as.matrix(W) %*% as.matrix(A)
	
	#RMSE <- Metrics::rmse(y_obs, y_hat)
	RMSE <- sqrt(mean((y_obs-y_hat)^2))
	
	return( RMSE )
}

ERA_tr <- function(y=y, X=X, nvar=nvar, dist=1, const = 0, lambda = 0, it = 0, ceps = 0.0001){
	
        if(dist==1) Yfam = "gaussian"
        
        data = cbind(y,X)
        ndset = length(nvar)
        ncase = nrow(y)
        ny = ncol(y)
        sum_nvar = sum(nvar)
        # ------------------
        # data normalization
        # ------------------
        dataN <- scale(data)/sqrt(ncase-1)
        y <- dataN[,1]
        X <- dataN[,-1]
        
        # ------------------
        # initialization
        # ------------------
        W0 <- matrix(0, nrow=sum_nvar, ncol=ndset)

        kk = 0
        for (j in 1:ndset){
                k = kk + 1
                kk = kk + nvar[,j]
                if ( nvar[,j]==1 ) {
                        W0[k:kk,j] = 1
                } else {
                        W0[k:kk,j] = 99
                }
        }
        windex <- which(W0 == 99)
        num_windex = length(windex)
        W = W0
        W[windex] <- runif(num_windex)
        F = X%*%W
        if( const==1 ) {
                demF <- chol(t(F)%*%F)
                F <- t(solve(t(demF),t(F)))
        } else {
                for (j in 1:ndset) { F[,j] <- F[,j] / norm(F[,j,drop=FALSE], type = "2") }
        }
        A <- solve(t(F)%*%F,t(F)%*%y)
        lp <- F%*%A
        adj.DV <- linkfun(X, lp, y, family = Yfam)
        mu <- adj.DV$mu #getmu(lp,dist)
        v <- adj.DV$mvar #getvmu(mu,dist)
        V <- diag(v)
        z <- adj.DV$z # lp + (y-mu)/v
        
        # ------------------
        # IRLS Algorithm
        # ------------------
        vecW <- W[windex]
        est_new <- c(vecW,A)
        est_old <- rep(0, length(est_new))

        while ( sum(abs(est_new - est_old)) > ceps ) {
                it = it + 1
                if ( it > 1000 ) {
                        print("Message: not converged within 1000 iterations")
                        break
                }
                est_old <- est_new
                
                # Step1: Update weights
                M0 <- kronecker(t(A),X)
                M <- M0[,windex]
                temp <- t(M)%*%V%*%M + lambda*diag(num_windex)
				try( vecW <- solve(temp, t(M)%*%V%*%z) )
				if(!is(vecW, "try-error")){
				vecW <- ginv(temp) %*% (t(M)%*%V%*%z)
				}
								
				W[windex] <- vecW
                F = X%*%W
                if( const==1 ) {
                        demF <- chol(t(F)%*%F)
                        F <- t(solve(t(demF),t(F)))
                } else {
                        for (j in 1:ndset) { F[,j] <- F[,j] / norm(F[,j,drop=FALSE], type = "2") }
                }
                
                # Step2: Update A
                A <- solve(t(F)%*%V%*%F, t(F)%*%V%*%z)
                
                # Step3: Update V and z
                lp <- F%*%A
                adj.DV <- linkfun(X, lp, y, family = Yfam)
                mu <- adj.DV$mu #getmu(lp,dist)
                v <- adj.DV$mvar #getvmu(mu,dist)
                V <- diag(v)
                z <- adj.DV$z # lp + (y-mu)/v
                est_new <- c(vecW,A)
        }
        vecW <- vecW
        A <- A

        output.ERA <- list(it=it, X=X, adj.DV=z, W=W, F=F, A=A)
        res <- output.ERA
        res$it
        res$X
        res$adj.DV
        res$W
        res$F
        res$A
        
        class(res) <- c("ERA_tr", class(res))
        invisible(res)
}

getdeviance <- function(y, mu, datatype) {
	if (datatype == 1) {DS = sum((y - mu)^2)} 
	return(DS)
}

linkfun <- function(x, lp, y, weights = NULL, offset = NULL, family = family,
					mustart = NULL, etastart = NULL, start = NULL)
{
	# Ref: stats::glm() in R
	# In R, family objects: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/glm.R
	# IWSL with adj DVs: https://www.statistics.ma.tum.de/fileadmin/w00bdb/www/czado/lec2.pdf
	# From linkfun(): adj.DV$w (weights in IWLS) is all 1 when gaussian
	
	x <- as.matrix(x)
	nvars <- ncol(x)
	nobs <- nrow(x)
	if (is.null(weights)) {	weights <- rep.int(1, nobs) }
    if (is.null(offset)) { offset <- rep.int(0, nobs) }
	
	## Define family:	
	if(is.character(family)) {
		family <- get(family, mode = "function", envir = parent.frame())
	}
	
	if(is.function(family)) family <- family()
	
	if(is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}
	
	## Get family functions:
	variance <- family$variance
	linkinv  <- family$linkinv
		 if (!is.function(variance) || !is.function(linkinv) ) {
			stop("'family' argument seems not to be a valid family object", call. = FALSE)
		}
    mu.eta <- family$mu.eta
	
	unless.null <- function(x, if.null) if(is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu  <- unless.null(family$validmu,  function(mu) TRUE)
	
	if(is.null(mustart)) {
        ## calculates mustart and may change y and weights and set n (!)
        eval(family$initialize)
    } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
	
	## Get (initial) eta
	coefold <- NULL
	eta <- lp
	mu <- linkinv(eta)
	mvar <- variance(mu) # the variance as a function of the mean
		
	if (!(validmu(mu) && valideta(eta))) {
		stop("cannot find valid starting values: please specify some", call. = FALSE)
	}

	## Calculate z and w
	good <- weights > 0
	mu.eta.val <- mu.eta(eta)
		if (any(is.na(mu.eta.val[good]))) {stop("NAs in d(mu)/d(eta)")}
	
	# drop observations for which w will be zero
	good <- (weights > 0) & (mu.eta.val != 0)
	z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
		# z: adjusted dependent variable
	w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
		# w: weights

	# BETA in GLM: Fisher Scoring
	# beta <- solve(t(X) %*% diag(w) %*% X) %*% (t(X) %*% diag(w) %*% z)
	
	# Outs
	return( list(mu=mu, mvar=mvar, z=z, w=w) )
}

