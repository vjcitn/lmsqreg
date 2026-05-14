
lmsqreg.fit <- function (YY, TT, edf = c(3, 5, 3), targlen = 50, targetx = seq(min(TT), 
    max(TT), length = targlen), pvec = c(0.05, 0.1, 0.25, 0.5, 
    0.75, 0.9, 0.95), maxit = 15, tol = 0.01, verb = FALSE, lam.fixed = NULL, 
    mu.fixed = NULL, sig.fixed = NULL, xcuts = quantile(TT, c(0.2, 
        0.4, 0.6, 0.8)), sig.init = mad(YY)/median(YY), lam.init = NULL) 
{
     yn <- deparse(substitute(YY))
     xn <- deparse(substitute(TT))
    qsys <- function(lmsobj, targetx, pvec = c(0.05, 0.1, 0.25, 
        0.5, 0.75, 0.9, 0.95)) {
        xrange <- range(lmsobj$ordt)
        Z <- qnorm(pvec)
        nro <- length(Z)
        outmat <- matrix(NA, nr = nro, nc = length(targetx))
        lmsmat <- cbind(sort(lmsobj$ordt), lmsobj$lam, lmsobj$mu, 
            lmsobj$sig)
        L <- approx(lmsmat[, 1], lmsmat[, 2], targetx, rule = 2)$y
        M <- approx(lmsmat[, 1], lmsmat[, 3], targetx, rule = 2)$y
        S <- approx(lmsmat[, 1], lmsmat[, 4], targetx, rule = 2)$y
        if (all(L != 0)) {
            for (i in 1:nro) outmat[i, ] <- M * (1 + L * S * 
                Z[i])^(1/L)
        }
        else if (all(L == 0)) {
            for (i in 1:nro) outmat[i, ] <- M * exp(S * Z[i])
        }
        dimnames(outmat) <- list(paste("P", as.character(pvec), 
            sep = ""), NULL)
        outlist <- list(outmat = outmat, targetx = targetx, pcts = pvec, 
            edf = lmsobj$edf)
        class(outlist) <- "qsys.out"
        outlist
    }
    validate <- function(qsys, t.val, y.val, rule = 2) {
        pout <- rep(NA, length(qsys$pcts))
        for (k in 1:length(qsys$pcts)) {
            pk <- approx(qsys$targetx, qsys$outmat[k, ], t.val, 
                rule = rule)
            pout[k] <- (sum(y.val < pk$y))/length(y.val)
        }
        list(pout, p.val = qsys$pcts)
    }
    fit.date <- date()
    fit.version <- Version(lmsqreg.fit)
    ot <- order(TT)
    TT <- TT[ot]
    YY <- YY[ot]
    N <- length(YY)
    lam <- rep(1, N)
    if (!is.null(lam.init)) {
        if (length(lam.init) == N) 
            lam <- lam.init
        else if (length(lam.init) == 1) 
            lam <- rep(lam.init, N)
        else {
            warning("lam.init not length 1 or N, using rep(lam.init[1],N)")
            lam <- rep(lam.init[1], N)
        }
    }
    if (!is.null(lam.fixed)) 
        lam <- rep(lam.fixed, N)
    Yshift <- 0
    if (any(YY < 1)) {
        Yshift <- 1 - min(YY)
        print(paste("Shifting Y by Yshift=", Yshift))
        YY <- YY + Yshift
    }
    mu <- fitted(loess(YY ~ TT))
    if (!is.null(mu.fixed)) 
        mu <- rep(mu.fixed, N)
    if (!is.null(sig.init)) {
        if (length(sig.init) == N) 
            sig <- sig.init
        else if (length(sig.init) == 1) 
            sig <- rep(sig.init, N)
        else {
            warning("sig.init not length 1 or N, using rep(sig.init[1],N)")
            sig <- rep(sig.init[1], N)
        }
    }
    else sig <- 1 + sqrt(pmax(lowess(pmax((YY/mu - 1), 0)^2)$y, 
        0))
    if (!is.null(sig.fixed)) 
        sig <- rep(sig.fixed, N)
    if (all(lam != 0)) 
        z <- ((YY/mu)^lam - 1)/(lam * sig)
    else if (all(lam == 0)) 
        z <- log(YY/mu)/sig
    else {
        z <- rep(NA, N)
        z[lam == 0] <- log(YY[lam == 0]/mu[lam == 0])/sig[lam == 
            0]
        z[lam != 0] <- ((YY[lam != 0]/mu[lam != 0])^lam[lam != 
            0] - 1)/(lam[lam != 0] * sig[lam != 0])
    }
    u <- function(y, lam, mu, sig) {
        YY <- y
        N <- length(y)
        if (all(lam != 0)) 
            z <- ((YY/mu)^lam - 1)/(lam * sig)
        else if (all(lam == 0)) 
            z <- log(YY/mu)/sig
        else {
            z <- rep(NA, N)
            z[lam == 0] <- log(YY[lam == 0]/mu[lam == 0])/sig[lam == 
                0]
            z[lam != 0] <- ((YY[lam != 0]/mu[lam != 0])^lam[lam != 
                0] - 1)/(lam[lam != 0] * sig[lam != 0])
        }
        z2m1 <- z * z - 1
        lyom <- log(y/mu)
        if (all(lam != 0)) 
            ul <- z/lam * (z - lyom/sig) - lyom * z2m1
        else if (all(lam == 0)) 
            ul <- rep(0, N)
        um <- z/(mu * sig) + (lam * z2m1)/mu
        us <- z2m1/sig
        list(lam = ul, mu = um, sig = us)
    }
    Wfun <- function(y, lam, mu, sig) {
        s2 <- sig * sig
        ls <- lam * sig
        l <- (7 * s2)/4
        m <- (1 + 2 * ls * ls)/(mu * mu * s2)
        s <- 2/s2
        lm <- -1/(2 * mu)
        ms <- (2 * lam)/(mu * sig)
        list(lam = l, mu = m, sig = s, lam.mu = lm, lam.sig = ls, 
            mu.sig = ms)
    }
    first <- TRUE
    iter <- 1
    Smooth.spline <- function(x, y, w, df) {
        #ans <- gam.spar.R(y ~ s(x, k = df, m = 2), w = w)
        #list(x = x, y = ans$fitted, rpen = ans$rpen)
        ss <- smooth.spline(x,y,w=w,df=df)
	ypred <- predict(ss,x)$y
        res <- y-ypred
        rpen <- t(ypred) %*% (w * res)
	list(x=x, y=ypred, rpen=rpen)
    }
    while ((first || nonconv) & iter < maxit) {
        iter <- iter + 1
        U <- u(YY, lam, mu, sig)
        W <- Wfun(YY, lam, mu, sig)
        if (first) {
            first <- FALSE
            psd1 <- U$lam/W$lam + lam
            nlam <- Smooth.spline(TT, psd1, w = W$lam, df = edf[1])$y
            psd2 <- U$mu/W$mu + mu - ((nlam - lam) * W$lam.mu)/W$mu
            nmu <- Smooth.spline(TT, psd2, w = W$mu, df = edf[2])$y
            psd3 <- U$sig/W$sig + sig - ((nlam - lam) * W$lam.sig)/W$sig - 
                ((nmu - mu) * W$mu.sig)/W$sig
            nsig <- Smooth.spline(TT, psd3, w = W$sig, df = edf[3])$y
        }
        psd1a <- U$lam/W$lam + nlam - ((nmu - mu) * W$lam.mu)/W$lam - 
            ((nsig - sig) * W$lam.sig)/W$lam
        lamtmp <- Smooth.spline(TT, psd1a, w = W$lam, df = edf[1])
        nlam <- lamtmp$y
        if (!is.null(lam.fixed)) 
            nlam <- rep(lam.fixed, N)
        psd2a <- U$mu/W$mu + mu - ((nlam - lam) * W$lam.mu)/W$mu - 
            ((nsig - sig) * W$lam.sig)/W$mu
        mutmp <- Smooth.spline(TT, psd2a, w = W$mu, df = edf[2])
        nmu <- mutmp$y
        if (!is.null(mu.fixed)) 
            nmu <- rep(mu.fixed, N)
        psd3a <- U$sig/W$sig + sig - ((nlam - lam) * W$lam.sig)/W$sig - 
            ((nmu - mu) * W$mu.sig)/W$sig
        sigtmp <- Smooth.spline(TT, psd3a, w = W$sig, df = edf[3])
        nsig <- sigtmp$y
        if (!is.null(sig.fixed)) 
            nsig <- rep(sig.fixed, N)
        if (verb) {
            print("lrange")
            print(range(lam - nlam))
            print("mrange")
            print(range(mu - nmu))
            print("srange")
            print(range(sig - nsig))
        }
        nonconv <- TRUE
        change <- max(c(abs(c(range(lam - nlam), range(mu - nmu), 
            range(sig - nsig)))))
        if (change < tol) 
            nonconv <- FALSE
        lam <- nlam
        mu <- nmu
        sig <- nsig
    }
    converged <- TRUE
    if (nonconv) {
        warning(paste(maxit, "iterations; did not converge, change=", 
            change))
        converged <- FALSE
    }
    if (all(lam != 0)) 
        z <- ((YY/mu)^lam - 1)/(lam * sig)
    else if (all(lam == 0)) 
        z <- log(YY/mu)/sig
    else {
        z <- rep(NA, N)
        z[lam == 0] <- log(YY[lam == 0]/mu[lam == 0])/sig[lam == 
            0]
        z[lam != 0] <- ((YY[lam != 0]/mu[lam != 0])^lam[lam != 
            0] - 1)/(lam[lam != 0] * sig[lam != 0])
    }
    rp.lam <- lamtmp$rpen
    if (!is.null(lam.fixed)) 
        rp.lam <- 0
    rp.mu <- mutmp$rpen
    if (!is.null(mu.fixed)) 
        rp.mu <- 0
    rp.sig <- sigtmp$rpen
    if (!is.null(sig.fixed)) 
        rp.sig <- 0
    upl <- sum(lam * log(YY/mu) - log(sig) - 0.5 * z * z)
    pl <- upl - 0.5 * rp.lam - 0.5 * rp.mu - 0.5 * rp.sig
    xfac <- cut(TT, round(c(min(TT) - 0.001, xcuts, max(TT) + 
        0.001), 3))
    zspl <- split(z, xfac)
    ntests <- length(zspl)
    ps <- rep(NA, ntests + 1)
    tps <- rep(NA, ntests + 1)
    vps <- rep(NA, ntests + 1)
    names(ps) <- c(names(table(xfac)), "Overall")
    names(tps) <- names(ps)
    names(vps) <- names(ps)
    unit.var.test <- function(x, nullv = 1) {
        V <- var(x)
        n <- length(x)
        if (V <= nullv) 
            ans <- (2 * pchisq(((n - 1) * V)/nullv, n - 1))
        else ans <- 2 * (1 - pchisq(((n - 1) * V)/nullv, n - 
            1))
        if (ans > 1) 
            ans <- 1
        list(var = V, p.val = ans)
    }
    for (i in 1:ntests) {
        ps[i] <- ks.test(zspl[[i]],"pnorm")$p.val
        tps[i] <- t.test(zspl[[i]])$p.val
        vps[i] <- unit.var.test(zspl[[i]])$p.val
    }
    ps[ntests + 1] <- ks.test(z,"pnorm")$p.val
    tps[ntests + 1] <- t.test(z)$p.val
    vps[ntests + 1] <- unit.var.test(z)$p.val
    lms.ans <- list(ordt = TT, lam = lam, mu = mu, sig = sig, 
        upl = upl, finalz = z, ps = ps, tps = tps, vps = vps, 
        edf = edf, niter = iter, converged = converged, fit.date = fit.date, 
        fitter.version = fit.version, yname = yn, xname = xn, 
        rpen = c(rp.lam, rp.mu, rp.sig), pl = pl, Yshift = Yshift)
    outq <- qsys(lms.ans, targetx, pvec)
    validout <- validate(outq, TT, YY)
    ans <- list(lms.ans = lms.ans, qsys = outq, validout = validout)
    class(ans) <- "lmsqreg.fit"
    ans
}

print.lmsqreg.fit <- function(x,...)
{
#Source version 2.5 97/03/26
#/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.print.lmsqreg.fit.q

	line1 <- paste("\nlms quantile regression, version ", x[[1]]$
		fitter.version, ", fit date ", x[[1]]$fit.date, "\n\n", sep = 
		"")
	cat(line1)
	line2 <- paste("Dependent variable:", x[[1]]$yname, 
		", independent variable:", x[[1]]$xname, "\n")
	cat(line2)
	if(x[[1]]$converged)
		line3 <- paste("The fit converged with EDF=(", paste(x[[1]]$edf,
			collapse = ","), "), PL=",round(x[[1]]$pl,3),"\n")
	else line3 <- paste("The fit failed to converge with EDF=(", paste(x[[1
			]]$edf, collapse = ","), "), after", x[[1]]$niter, 
			"iterations.\n")
	cat(line3)
	vout <- x[[3]]
	nom <- vout$p.val
	obs <- vout[[1]]
	val <- rbind(nom, obs)
	dimnames(val) <- list(c("nominal percentile", "estimated percentile"), 
		rep(" ", length(nom)))
	print(round(val, 3))
	cat("\nKS tests: (intervals in",x[[1]]$xname,"//p-values)\n")
	print(round(x[[1]]$ps,3))
	cat("\nt tests: (intervals in",x[[1]]$xname,"//p-values)\n")
	print(round(x[[1]]$tps,3))
	cat("\nX2 tests (unit variance): (intervals in",x[[1]]$xname,"//p-values)\n")
	print(round(x[[1]]$vps,3))
	invisible(0)
}
plot.lmsqreg.fit <- function (x, fullsys = TRUE, medname = "P0.5", xlab = " ", ylab = " ", 
    Yshift = -x[[1]]$Yshift, title = paste("LMS fit with edf = (", 
        ob$edf[1], ",", ob$edf[2], ",", ob$edf[3], "), PL=", 
        round(x[[1]]$pl, 3), sep = ""), CEX = 1.1, tx = function(z) z, 
    xaxat = NULL, xaxlab = NULL, ...) 
{
    nnf <- matrix(c(1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4), nr = 4, 
        byrow = TRUE)
    layout(nnf, heights = c(1.5, 1, 1, 1))
    obj <- x
    ob <- obj$qsys
    YL <- range(ob$outmat) + Yshift
    XL <- range(ob$targetx)
    YLAB <- obj[[1]]$yname
    XLAB <- obj[[1]]$xname
    plot(tx(obj[[1]]$ordt), obj[[1]]$lam, xlab = obj[[1]]$xname, 
        ylab = "Lambda", cex = CEX, axes = FALSE)
    axis(2)
    if (is.null(xaxat)) 
        axis(1)
    else axis(1, at = xaxat, lab = xaxlab)
    plot(tx(obj[[1]]$ordt), obj[[1]]$mu + Yshift, xlab = obj[[1]]$xname, 
        ylab = "Mu", cex = CEX, axes = FALSE)
    axis(2)
    if (is.null(xaxat)) 
        axis(1)
    else axis(1, at = xaxat, lab = xaxlab)
    plot(tx(obj[[1]]$ordt), obj[[1]]$sig, xlab = obj[[1]]$xname, 
        ylab = "Sigma", cex = CEX, axes = FALSE)
    axis(2)
    if (is.null(xaxat)) 
        axis(1)
    else axis(1, at = xaxat, lab = xaxlab)
    if (ylab != " ") 
        YLAB <- ylab
    if (xlab != " ") 
        XLAB <- xlab
    plot(tx(ob$targetx), ob$outmat[medname, ] + Yshift, type = "l", 
        xlab = XLAB, ylab = YLAB, xlim = tx(XL), ylim = YL, main = title, 
        cex = CEX, axes = FALSE)
    axis(2)
    if (is.null(xaxat)) 
        axis(1)
    else axis(1, at = xaxat, lab = xaxlab)
    for (j in 1:nrow(ob$outmat)) {
        text(tx(XL[2]), ob$outmat[j, ncol(ob$outmat)] + Yshift, 
            as.character(ob$pcts[j]), cex = CEX)
        lines(tx(ob$targetx), ob$outmat[j, ] + Yshift, lty = 2, 
            cex = CEX)
    }
    invisible(0)
}
zscores <- function(y, x, obj)
{
# create zscores for arbitrary (y,x) given quantile regression object obj
#Source version 1.4 97/03/26
#/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.zscores.q
# revision prior to camacs tutorial
	if(class(obj) != "lmsqreg.fit") stop("obj must be class lmsqreg.fit")
	baserange <- range(obj[[1]]$ordt)
	topin <- max(x)
	botin <- min(x)
	if(topin > baserange[2])
		warning(paste("constant extrap. from", round(baserange[2], 4), 
			"to", round(topin, 4)))
	if(botin < baserange[1]) warning(paste("constant extrap. from", round(
			baserange[1], 4), "to", round(botin, 4)))	
	# interpolate fitted functions to target x values
	lapp <- approx(obj[[1]]$ordt, obj[[1]]$lam, x, rule = 2)$y
	mapp <- approx(obj[[1]]$ordt, obj[[1]]$mu, x, rule = 2)$y
	sapp <- approx(obj[[1]]$ordt, obj[[1]]$sig, x, rule = 2)$y
	if(all(obj[[1]]$lam != 0))
		z <- ((y/mapp)^lapp - 1)/(lapp * sapp)
	else if(all(obj[[1]]$lam == 0))
		z <- log(y/mapp)/sapp
	else {
		bad <- obj[[1]]$lam == 0
		z <- lapp - lapp
		z[bad] <- log(y[bad]/mapp[bad])/sapp[bad]
		z[ - bad] <- ((y[ - bad]/mapp[ - bad])^lapp[ - bad] - 1)/(lapp[ -
			bad] * sapp[ - bad])
	}
	z
}
local.winsorization <- function(x, y, ncut = 5, k = 20)
{
require(parody)
# SCCS version local.winsorization.q 1.1 97/04/24
# /usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.local.winsorization.q
# this function winsorizes outliers
# in a step-function regression
#	lamtab <- function(n, k, alpha = 0.05)
#	{
## obtain the sequence of k critical values
## to be used in checking for k outliers
#		l <- 0:(k - 1)
#		p <- 1 - ((alpha/2)/(n - l))	
#	# note, this can trip a bug in Splus 3.2 qt if n is very large (>50000?)
## the invibeta routine needs to be iterated some more
#		tvals <- qt(df = n - l - 2, p = p)
#		lams <- ((n - l - 1) * tvals)/sqrt((n - l - 2 + tvals^2) * (n - 
#			l))
#		lams
#	}
#	assign("lamtab", lamtab, f = 0)
#	skesd <- function(x, k)
#	{
## obtain k large abs. studentized residuals
#		bigres <- rep(NA, k)
#		bigresind <- rep(NA, k)
#		n <- length(x)
#		curx <- x
#		ind <- 1:n
#		inds <- NULL
#		for(i in 1:k) {
#			last <- n - i + 1
#			m <- mean(curx)
#			ares <- abs(curx - m)
#			oares <- order(ares)
#			bigres[i] <- ares[oares[last]]/sqrt(var(curx))
#			curx <- curx[ - oares[last]]
#			inds <- c(inds, ind[oares[last]])
#			ind <- ind[ - oares[last]]
#		}
#		if(length(inds) == 0) {
#			inds <- NA
#			res <- NA
#		}
#		list(res = bigres, ind = inds)
#	}
#	assign("skesd", skesd, f = 0)
#	sgesdri <- function(x, k = ((length(x) %% 2) * floor(length(x)/2) + (1 - (
#		length(x) %% 2)) * (length(x)/2 - 1)), alpha = 0.05)
#	{
## obtain Rosner GESD outliers and indices
#		E <- skesd(x, k)
#		R <- E$res
#		I <- E$ind
#		n <- length(x)	#if(n == 100 & k == 49)
#		L <- lamtab(length(x), k, alpha = alpha)
#		worst <- max((1:k)[R > L])
#		if(is.na(worst) | is.null(worst))
#			list(ind = NA, val = NA)
#		else list(ind = I[1:worst], val = x[I[1:worst]])
#	}
#	assign("sgesdri", sgesdri, f = 0)
# begin the local winsorization -- obtain a crude
# cut data as requested and seek outliers in each section
	cutfs <- (0:ncut)/ncut
	lims <- quantile(x, cutfs)
	lims[1] <- lims[1] - 0.1
	lims[length(lims)] <- lims[length(lims)] + 0.1
	xc <- cut(x, lims)	#
	iinds <- 1:length(x)
	iindsspl <- split(iinds, xc)
	xspl <- split(x, xc)
	yspl <- split(y, xc)
	nout <- list()
	bad <- list()
	for(i in 1:length(yspl)) {
		curo <- calout.detect(yspl[[i]], meth="GESD")$ind
		if(length(curo) >= 1 & !is.na(curo[1])) {
			bad[[i]] <- iindsspl[[i]][curo]
			nout[[i]] <- length(curo)
			curcln <- yspl[[i]][ - curo]
			lowcln <- min(curcln)
			hicln <- max(curcln)
			curbad <- yspl[[i]][curo]
			for(j in 1:nout[[i]]) {
				if(curbad[j] > mean(curcln))
				  yspl[[i]][curo[j]] <- hicln
				else yspl[[i]][curo[j]] <- lowcln
			}
		}
	}
	if(sum(unlist(nout)) == 0) {
		print("no local outliers")
		return(list(x = x, y = y))	# no outliers;
	}
	else print(paste(sum(unlist(nout)), "outliers"))
	retord <- order(as.integer(unlist(iindsspl)))
	return(list(x = as.double(unlist(xspl))[retord], y = as.double(unlist(yspl))[
		retord], bad = unlist(bad)))
}
get.quantiles <- function(lmsobj, targetx, pvec = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
{
# qsys: function to develop the quantile estimates
# from an LMS fit
# get.quantiles.q 1.1 97/04/24
# /usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.get.quantiles.q
	if(class(lmsobj) != "lmsqreg.fit") stop(
			"only applies to lmsqreg.fit object")
	lmsobj <- lmsobj[[1]]
	xrange <- range(lmsobj$ordt)
	Z <- qnorm(pvec)
	nro <- length(Z)
	outmat <- matrix(NA, nr = nro, nc = length(targetx))
	lmsmat <- cbind(sort(lmsobj$ordt), lmsobj$lam, lmsobj$mu, lmsobj$sig)
	L <- approx(lmsmat[, 1], lmsmat[, 2], targetx, rule = 2)$y
	M <- approx(lmsmat[, 1], lmsmat[, 3], targetx, rule = 2)$y
	S <- approx(lmsmat[, 1], lmsmat[, 4], targetx, rule = 2)$y
	if(all(L != 0)) {
		for(i in 1:nro)
			outmat[i,  ] <- M * (1 + L * S * Z[i])^(1/L)
	}
	else if(all(L == 0)) {
		for(i in 1:nro)
			outmat[i,  ] <- M * exp(S * Z[i])
	}
	dimnames(outmat) <- list(paste("P", as.character(pvec), sep = ""), NULL
		)
	outlist <- list(outmat = outmat, targetx = targetx, pcts = pvec, edf = 
		lmsobj$edf)
	class(outlist) <- "qsys.out"
	outlist
}
#gam.spar.R <- function (formula, weights = NULL) 
#{
#    #require(mgcv)
#    #fit <- gam(formula)
#    require(modreg)
#    w <- weights
#    if (!length(w)) 
#        w <- rep(1, nrow(m))
#    else if (any(w < 0)) 
#        stop("negative weights not allowed")
#    rpen <- t(fit$fitted) %*% (w * fit$residuals)
#    list(y = fit$fitted, rpen = rpen, fitted = fit$fitted)
#}
lmsqreg.search <- function(y, x, startedf = c(3, 5, 3), boundedf = c(1, 2, 1), search.seq = c(1,3,2),...)
{
# source 1.1 97/05/25
# /usr16/stdevs/stdev0f/SLIBS/lmsqregdev/SCCS/s.lmsqreg.search.q
	oldfit <- lmsqreg.fit(y, x, edf = startedf, ...)
	oldpl <- oldfit[[1]]$pl
	print(oldpl)
	newedf <- curedf <- startedf
	newpl <- oldpl
	for(edfcomp in search.seq) {
		while(1) {
			newedf[edfcomp] <- curedf[edfcomp] - 1
			if(newedf[edfcomp] < boundedf[edfcomp]) {
				newedf <- curedf
				break
			}
			print(newedf)
			tmpfit <- lmsqreg.fit(y, x, edf = newedf, ...)
			newpl <- tmpfit[[1]]$pl
			print(newpl)
			if(newpl < (oldpl - 2))
				break
			else {
				oldfit <- tmpfit
				oldpl <- newpl
				curedf <- newedf
			}
		}
	}
	list(fit = oldfit, pl = oldpl, edf = curedf)
}
locwin.envelope <- function(x, y, ncut = 5, k = 20)
{
# this function displays winsorization bands
# support functions
	lamtab <- function(n, k, alpha = 0.05)
	{
# obtain the sequence of k critical values
# to be used in checking for k outliers
		l <- 0:(k - 1)
		p <- 1 - ((alpha/2)/(n - l))	
	# note, this can trip a bug in Splus 3.2 qt if n is very large (>50000?)
# the invibeta routine needs to be iterated some more
		tvals <- qt(df = n - l - 2, p = p)
		lams <- ((n - l - 1) * tvals)/sqrt((n - l - 2 + tvals^2) * (n - 
			l))
		lams
	}
	assign("lamtab", lamtab, f = 0)
	skesd <- function(x, k)
	{
# obtain k large abs. studentized residuals
		bigres <- rep(NA, k)
		bigresind <- rep(NA, k)
		n <- length(x)
		curx <- x
		ind <- 1:n
		inds <- NULL
		for(i in 1:k) {
			last <- n - i + 1
			m <- mean(curx)
			ares <- abs(curx - m)
			oares <- order(ares)
			bigres[i] <- ares[oares[last]]/sqrt(var(curx))
			curx <- curx[ - oares[last]]
			inds <- c(inds, ind[oares[last]])
			ind <- ind[ - oares[last]]
		}
		if(length(inds) == 0) {
			inds <- NA
			res <- NA
		}
		list(res = bigres, ind = inds)
	}
	sgesdri <- function(x, k = ((length(x) %% 2) * floor(length(x)/2) + (1 - (
		length(x) %% 2)) * (length(x)/2 - 1)), alpha = 0.05)
	{
# obtain Rosner GESD outliers and indices
		E <- skesd(x, k)
		R <- E$res
		I <- E$ind
		n <- length(x)	#if(n == 100 & k == 49)
		L <- lamtab(length(x), k, alpha = alpha)
		worst <- max((1:k)[R > L])
		if(is.na(worst) | is.null(worst))
			list(ind = NA, val = NA)
		else list(ind = I[1:worst], val = x[I[1:worst]])
	}
# begin the local winsorization -- obtain a crude
# cut data as requested and seek outliers in each section
	cutfs <- (0:ncut)/ncut
	lims <- quantile(x, cutfs)
	lims[1] <- lims[1] - 0.1
	lims[length(lims)] <- lims[length(lims)] + 0.1
	xc <- cut(x, lims)	#
	iinds <- 1:length(x)
	iindsspl <- split(iinds, xc)
	xspl <- split(x, xc)
	yspl <- split(y, xc)
	yspl2 <- list()
	nout <- list()
	bad <- list()
	for(i in 1:length(yspl)) {
		curo <- sgesdri(yspl[[i]])$ind
		bad[[i]] <- iindsspl[[i]][curo]
		if(!is.na(curo[1])) {
			nout[[i]] <- length(curo)
			curcln <- yspl[[i]][ - curo]
		}
		else {
			nout[[i]] <- 0
			curcln <- yspl[[i]]
		}
		lowcln <- min(curcln)
		hicln <- max(curcln)
		yspl[[i]] <- rep(hicln, length(yspl[[i]]))
		yspl2[[i]] <- rep(lowcln, length(yspl[[i]]))
	}
	print(paste(sum(unlist(nout)), "outliers"))
	retord <- order(as.integer(unlist(iindsspl)))
	return(x = as.double(unlist(xspl))[retord], y = as.double(unlist(yspl))[
		retord], y2 = as.double(unlist(yspl2)[retord]), bad = unlist(
		bad))
}
validate.report <- function(qreg, y.val, t.val, 
        rule = 2, xcuts = quantile(t.val, c(0.2, 0.4, 0.6, 
	0.8)))
{
# validate.report.q 1.1 97/04/24
# /usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.validate.report.q
# function to validate a quantile regression system
# on new data at t.val, y.val
	qsys <- qreg$qsys
	pout <- rep(NA, length(qsys$pcts))
	for(k in 1:length(qsys$pcts)) {
		pk <- approx(qsys$targetx, qsys$outmat[k,  ], t.val, rule = 
			rule)
		pout[k] <- (sum(y.val < pk$y))/length(y.val)
	}
	z <- zscores(y.val, t.val, qreg)
	TT <- t.val
	xcuts <- quantile(t.val, c(0.2, 0.4, 0.6, 0.8))
	xfac <- cut(TT, round(c(min(TT) - 0.001, xcuts, max(TT) + 0.001), 3))
	zspl <- split(z, xfac)
	ntests <- length(zspl)
	ps <- rep(NA, ntests + 1)
	tps <- rep(NA, ntests + 1)
	vps <- rep(NA, ntests + 1)
	names(ps) <- c(names(table(xfac)), "Overall")
	names(tps) <- names(ps)
	names(vps) <- names(ps)
	unit.var.test <- function(x, nullv = 1)
	{
		V <- var(x)
		n <- length(x)
		if(V <= nullv)
			ans <- (2 * pchisq(((n - 1) * V)/nullv, n - 1))
		else ans <- 2 * (1 - pchisq(((n - 1) * V)/nullv, n - 1))
		if (ans > 1) ans <- 1
		list(var = V, p.val = ans)
	}
	for(i in 1:ntests) {
		ps[i] <- ks.test(zspl[[i]],"pnorm")$p.val
		tps[i] <- t.test(zspl[[i]])$p.val
		vps[i] <- unit.var.test(zspl[[i]])$p.val
	}
	ps[ntests + 1] <- ks.test(z,"pnorm")$p.val
	tps[ntests + 1] <- t.test(z)$p.val
	vps[ntests + 1] <- unit.var.test(z)$p.val
	cat("Nominal quantile coverage\n")
	cat(round(qsys$pcts, 5))
	cat("\n")
	cat("Actual quantile coverage\n")
	cat(round(pout, 5))
	cat("\n")
	cat("\nKS tests: (intervals in", qreg[[1]]$xname, 
		"//p-values)\n")
	print(round(ps, 3))
	cat("\nt tests: (intervals in", qreg[[1]]$xname, "//p-values)\n")
	print(round(tps, 3))
	cat("\nX2 tests (unit variance): (intervals in", qreg[[1]]$xname, 
		"//p-values)\n")
	print(round(vps, 3))
	invisible(list(pout, p.val = qsys$pcts, ps = ps, tps = tps, vps = vps))
}
Version <- function(x)
attr(x, "version")


