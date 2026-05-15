## Suppress R CMD CHECK notes for ggplot2 aesthetics column names
utils::globalVariables(c("x", "y", "pct", "label"))

#' Fit LMS quantile regression
#'
#' Fits the LMS (Lambda-Mu-Sigma) quantile regression model of Cole and Green
#' (1992) using penalized likelihood with cubic smoothing splines.
#'
#' @param YY numeric response vector (must be positive; shifted automatically if
#'   any values are below 1)
#' @param TT numeric predictor vector
#' @param edf length-3 numeric vector of equivalent degrees of freedom for the
#'   lambda, mu, and sigma smoothing splines
#' @param targlen integer; number of points in the target x grid
#' @param targetx numeric vector of x values at which quantiles are evaluated
#' @param pvec numeric vector of probabilities for quantile output
#' @param maxit integer; maximum number of Fisher-scoring iterations
#' @param tol numeric; convergence tolerance (max absolute parameter change)
#' @param verb logical; if `TRUE`, print parameter-change ranges each iteration
#' @param lam.fixed numeric scalar; if not `NULL`, hold lambda fixed at this
#'   value throughout
#' @param mu.fixed numeric scalar; if not `NULL`, hold mu fixed at this value
#' @param sig.fixed numeric scalar; if not `NULL`, hold sigma fixed at this
#'   value
#' @param xcuts numeric vector of x-quantile cut points used to define
#'   intervals for goodness-of-fit tests
#' @param sig.init initial sigma value(s); length 1 or `length(YY)`
#' @param lam.init initial lambda value(s); length 1 or `length(YY)`
#'
#' @return An object of class `"lmsqreg.fit"`, a list with three components:
#'   \describe{
#'     \item{lms.ans}{list of fitted LMS curves and convergence diagnostics}
#'     \item{qsys}{quantile system, class `"qsys.out"`, with `outmat`,
#'       `targetx`, `pcts`, and `edf`}
#'     \item{validout}{in-sample quantile coverage from `validate()`}
#'   }
#'
#' @references Cole, T.J. and Green, P.J. (1992). Smoothing reference centile
#'   curves: The LMS method and penalized likelihood. *Statistics in Medicine*,
#'   **11**, 1305–1319.
#'
#' @examples
#' set.seed(123)
#' x <- runif(300, 10, 20)
#' y <- 8 + 2 * sin(x) + rnorm(300, 0, x / 11)
#' fit <- lmsqreg.fit(y, x, maxit = 25)
#' print(fit)
#'
#' @importFrom stats approx smooth.spline predict fitted loess lowess median
#' @importFrom stats ks.test t.test pchisq qnorm qt quantile var mad
#' @export
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
        outmat <- matrix(NA, nrow = nro, ncol = length(targetx))
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
        message(paste("Shifting Y by Yshift=", Yshift))
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
        ss <- smooth.spline(x, y, w = w, df = df)
        ypred <- predict(ss, x)$y
        res <- y - ypred
        rpen <- t(ypred) %*% (w * res)
        list(x = x, y = ypred, rpen = rpen)
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
            message("lrange: ", paste(range(lam - nlam), collapse = " "))
            message("mrange: ", paste(range(mu - nmu), collapse = " "))
            message("srange: ", paste(range(sig - nsig), collapse = " "))
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
            change, "; consider increasing maxit"))
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
        ps[i] <- ks.test(zspl[[i]], "pnorm")$p.val
        tps[i] <- t.test(zspl[[i]])$p.val
        vps[i] <- unit.var.test(zspl[[i]])$p.val
    }
    ps[ntests + 1] <- ks.test(z, "pnorm")$p.val
    tps[ntests + 1] <- t.test(z)$p.val
    vps[ntests + 1] <- unit.var.test(z)$p.val
    lms.ans <- list(ordt = TT, rawY = YY, lam = lam, mu = mu, sig = sig,
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

#' Print an lmsqreg.fit object
#'
#' Displays convergence status, equivalent degrees of freedom, penalized
#' log-likelihood, nominal vs. estimated quantile coverage, and
#' goodness-of-fit test p-values (KS, t, and unit-variance chi-squared)
#' within x-intervals.
#'
#' @param x an object of class `"lmsqreg.fit"`
#' @param ... further arguments (currently unused)
#' @return Invisibly returns `0`.
#' @method print lmsqreg.fit
#' @export
print.lmsqreg.fit <- function(x, ...)
{
    line1 <- paste("\nlms quantile regression, version ", x[[1]]$
        fitter.version, ", fit date ", x[[1]]$fit.date, "\n\n", sep =
        "")
    cat(line1)
    line2 <- paste("Dependent variable:", x[[1]]$yname,
        ", independent variable:", x[[1]]$xname, "\n")
    cat(line2)
    if (x[[1]]$converged)
        line3 <- paste("The fit converged with EDF=(", paste(x[[1]]$edf,
            collapse = ","), "), PL=", round(x[[1]]$pl, 3), "\n")
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
    cat("\nKS tests: (intervals in", x[[1]]$xname, "//p-values)\n")
    print(round(x[[1]]$ps, 3))
    cat("\nt tests: (intervals in", x[[1]]$xname, "//p-values)\n")
    print(round(x[[1]]$tps, 3))
    cat("\nX2 tests (unit variance): (intervals in", x[[1]]$xname, "//p-values)\n")
    print(round(x[[1]]$vps, 3))
    invisible(0)
}

#' Plot an lmsqreg.fit object
#'
#' Produces a four-panel ggplot2 display: the fitted lambda, mu, and sigma
#' curves (top row) and the full centile curve system (bottom, larger panel).
#'
#' @param x an object of class `"lmsqreg.fit"`
#' @param fullsys logical; reserved for future use
#' @param medname character; row name of the median in the quantile matrix
#' @param xlab character; x-axis label (defaults to the predictor variable
#'   name stored in `x`)
#' @param ylab character; y-axis label (defaults to the response variable name)
#' @param Yshift numeric; shift added back to y for display (default reverses
#'   any automatic shift applied during fitting)
#' @param title character; main title for the centile panel
#' @param CEX numeric; size scaling factor for points and text
#' @param show.data logical; if `TRUE` (default), overlay the raw data points
#'   in grey on the centile panel
#' @param tx function applied to x values before plotting (e.g. for a date
#'   transformation)
#' @param xaxat numeric vector of custom x-axis tick positions
#' @param xaxlab character vector of custom x-axis tick labels (parallel to
#'   `xaxat`)
#' @param ... further arguments (currently unused)
#' @return A `patchwork` ggplot2 object, printed for its side effect.
#' @method plot lmsqreg.fit
#' @import ggplot2
#' @import patchwork
#' @export
plot.lmsqreg.fit <- function(x, fullsys = TRUE, medname = "P0.5",
    xlab = NULL, ylab = NULL,
    Yshift = -x[[1]]$Yshift,
    title = NULL, show.data = TRUE,
    CEX = 1.1, tx = function(z) z,
    xaxat = NULL, xaxlab = NULL, ...)
{
    ob <- x$qsys
    obj <- x

    YLAB <- if (!is.null(ylab)) ylab else obj[[1]]$yname
    XLAB <- if (!is.null(xlab)) xlab else obj[[1]]$xname

    if (is.null(title)) {
        title <- paste0("LMS fit with edf = (", ob$edf[1], ",",
            ob$edf[2], ",", ob$edf[3], "), PL=",
            round(x[[1]]$pl, 3))
    }

    ordt_tx   <- tx(obj[[1]]$ordt)
    targetx_tx <- tx(ob$targetx)

    x_scale <- if (!is.null(xaxat))
        scale_x_continuous(breaks = xaxat, labels = xaxlab)
    else
        scale_x_continuous()

    p_lam <- ggplot(
        data.frame(x = ordt_tx, y = obj[[1]]$lam),
        aes(x = x, y = y)
    ) +
        geom_point(size = CEX * 0.5) +
        labs(x = obj[[1]]$xname, y = "Lambda") +
        x_scale +
        theme_bw()

    p_mu <- ggplot(
        data.frame(x = ordt_tx, y = obj[[1]]$mu + Yshift),
        aes(x = x, y = y)
    ) +
        geom_point(size = CEX * 0.5) +
        labs(x = obj[[1]]$xname, y = "Mu") +
        x_scale +
        theme_bw()

    p_sig <- ggplot(
        data.frame(x = ordt_tx, y = obj[[1]]$sig),
        aes(x = x, y = y)
    ) +
        geom_point(size = CEX * 0.5) +
        labs(x = obj[[1]]$xname, y = "Sigma") +
        x_scale +
        theme_bw()

    nq <- nrow(ob$outmat)
    centile_df <- do.call(rbind, lapply(seq_len(nq), function(j) {
        data.frame(
            x     = targetx_tx,
            y     = ob$outmat[j, ] + Yshift,
            pct   = as.character(ob$pcts[j]),
            stringsAsFactors = FALSE
        )
    }))
    label_df <- data.frame(
        x     = max(targetx_tx),
        y     = ob$outmat[, ncol(ob$outmat)] + Yshift,
        label = as.character(ob$pcts),
        stringsAsFactors = FALSE
    )

    data_df <- data.frame(
        x = ordt_tx,
        y = obj[[1]]$rawY + Yshift
    )
    p_centile <- ggplot(centile_df, aes(x = x, y = y, group = pct)) +
        geom_line(linetype = "dashed") +
        geom_text(
            data = label_df,
            aes(x = x, y = y, label = label),
            hjust = -0.1, size = CEX * 2.5,
            inherit.aes = FALSE
        ) +
        labs(x = XLAB, y = YLAB, title = title) +
        x_scale +
        theme_bw() +
        theme(legend.position = "none")
    if (show.data)
        p_centile <- p_centile +
            geom_point(data = data_df, aes(x = x, y = y),
                colour = "grey60", size = CEX * 0.4, inherit.aes = FALSE)

    (p_lam | p_mu | p_sig) / p_centile +
        plot_layout(heights = c(1.5, 3))
}

#' Compute LMS z-scores for new data
#'
#' Transforms arbitrary `(y, x)` pairs to standard-normal z-scores using the
#' LMS functions from a fitted `lmsqreg.fit` object. Linear interpolation is
#' used within the fitted x range; constant extrapolation is applied outside
#' it (with a warning).
#'
#' @param y numeric vector of response values
#' @param x numeric vector of predictor values (same length as `y`)
#' @param obj an object of class `"lmsqreg.fit"`
#' @return numeric vector of z-scores of the same length as `y`
#' @examples
#' set.seed(123)
#' xf <- runif(200, 10, 20)
#' yf <- 8 + 2 * sin(xf) + rnorm(200, 0, xf / 11)
#' fit <- lmsqreg.fit(yf, xf)
#' xnew <- runif(50, 10, 20)
#' ynew <- 8 + 2 * sin(xnew) + rnorm(50, 0, xnew / 11)
#' z <- zscores(ynew, xnew, fit)
#' cat("z-score mean:", round(mean(z), 3), "\n")
#' @export
zscores <- function(y, x, obj)
{
    if (!inherits(obj, "lmsqreg.fit")) stop("obj must be class lmsqreg.fit")
    baserange <- range(obj[[1]]$ordt)
    topin <- max(x)
    botin <- min(x)
    if (topin > baserange[2])
        warning(paste("constant extrap. from", round(baserange[2], 4),
            "to", round(topin, 4)))
    if (botin < baserange[1])
        warning(paste("constant extrap. from", round(baserange[1], 4),
            "to", round(botin, 4)))
    lapp <- approx(obj[[1]]$ordt, obj[[1]]$lam, x, rule = 2)$y
    mapp <- approx(obj[[1]]$ordt, obj[[1]]$mu,  x, rule = 2)$y
    sapp <- approx(obj[[1]]$ordt, obj[[1]]$sig, x, rule = 2)$y
    if (all(obj[[1]]$lam != 0))
        z <- ((y/mapp)^lapp - 1)/(lapp * sapp)
    else if (all(obj[[1]]$lam == 0))
        z <- log(y/mapp)/sapp
    else {
        bad <- obj[[1]]$lam == 0
        z <- lapp - lapp
        z[bad]  <- log(y[bad]/mapp[bad])/sapp[bad]
        z[-bad] <- ((y[-bad]/mapp[-bad])^lapp[-bad] - 1)/(lapp[-bad] * sapp[-bad])
    }
    z
}

#' Local winsorization of outliers
#'
#' Divides the data into `ncut` quantile-based sections and applies the GESD
#' (Generalized Extreme Studentized Deviate) outlier test within each section.
#' Detected outliers are replaced ("winsorized") by the nearest non-outlier
#' extreme of the same section.
#'
#' @param x numeric predictor vector
#' @param y numeric response vector (same length as `x`)
#' @param ncut integer; number of quantile-based sections
#' @param k integer; maximum number of outliers tested per section
#' @return A list with components:
#'   \describe{
#'     \item{x}{predictor values (unchanged)}
#'     \item{y}{winsorized response values}
#'     \item{bad}{integer indices of replaced observations (only present when
#'       outliers were found)}
#'   }
#' @references Rosner, B. (1983). Percentage points for a generalized ESD
#'   many-outlier procedure. *Technometrics*, **25**, 165–172.
#' @importFrom parody calout.detect
#' @export
local.winsorization <- function(x, y, ncut = 5, k = 20)
{
    cutfs <- (0:ncut)/ncut
    lims <- quantile(x, cutfs)
    lims[1] <- lims[1] - 0.1
    lims[length(lims)] <- lims[length(lims)] + 0.1
    xc <- cut(x, lims)
    iinds <- 1:length(x)
    iindsspl <- split(iinds, xc)
    xspl <- split(x, xc)
    yspl <- split(y, xc)
    nout <- list()
    bad  <- list()
    for (i in 1:length(yspl)) {
        curo <- calout.detect(yspl[[i]], method = "GESD")$ind
        if (length(curo) >= 1 & !is.na(curo[1])) {
            bad[[i]]  <- iindsspl[[i]][curo]
            nout[[i]] <- length(curo)
            curcln  <- yspl[[i]][-curo]
            lowcln  <- min(curcln)
            hicln   <- max(curcln)
            curbad  <- yspl[[i]][curo]
            for (j in 1:nout[[i]]) {
                if (curbad[j] > mean(curcln))
                    yspl[[i]][curo[j]] <- hicln
                else
                    yspl[[i]][curo[j]] <- lowcln
            }
        }
    }
    if (sum(unlist(nout)) == 0) {
        message("no local outliers")
        return(list(x = x, y = y))
    }
    message(paste(sum(unlist(nout)), "outliers"))
    retord <- order(as.integer(unlist(iindsspl)))
    return(list(
        x   = as.double(unlist(xspl))[retord],
        y   = as.double(unlist(yspl))[retord],
        bad = unlist(bad)
    ))
}

#' Extract quantile estimates from an LMS fit
#'
#' Replicates the internal `qsys()` computation outside of `lmsqreg.fit()`,
#' allowing quantiles to be evaluated at a user-supplied grid of x values.
#'
#' @param lmsobj an object of class `"lmsqreg.fit"`
#' @param targetx numeric vector of x values at which to evaluate quantiles
#' @param pvec numeric vector of probabilities
#' @return An object of class `"qsys.out"` with components `outmat`
#'   (matrix of quantile values, rows = quantiles, cols = `targetx`),
#'   `targetx`, `pcts`, and `edf`.
#' @export
get.quantiles <- function(lmsobj, targetx,
    pvec = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
{
    if (!inherits(lmsobj, "lmsqreg.fit"))
        stop("only applies to lmsqreg.fit object")
    lmsobj <- lmsobj[[1]]
    xrange <- range(lmsobj$ordt)
    Z <- qnorm(pvec)
    nro <- length(Z)
    outmat <- matrix(NA, nrow = nro, ncol = length(targetx))
    lmsmat <- cbind(sort(lmsobj$ordt), lmsobj$lam, lmsobj$mu, lmsobj$sig)
    L <- approx(lmsmat[, 1], lmsmat[, 2], targetx, rule = 2)$y
    M <- approx(lmsmat[, 1], lmsmat[, 3], targetx, rule = 2)$y
    S <- approx(lmsmat[, 1], lmsmat[, 4], targetx, rule = 2)$y
    if (all(L != 0)) {
        for (i in 1:nro)
            outmat[i, ] <- M * (1 + L * S * Z[i])^(1/L)
    }
    else if (all(L == 0)) {
        for (i in 1:nro)
            outmat[i, ] <- M * exp(S * Z[i])
    }
    dimnames(outmat) <- list(paste("P", as.character(pvec), sep = ""), NULL)
    outlist <- list(outmat = outmat, targetx = targetx, pcts = pvec,
        edf = lmsobj$edf)
    class(outlist) <- "qsys.out"
    outlist
}

#' Grid search over EDF parameters
#'
#' Optimizes the equivalent degrees of freedom (EDF) triple for an LMS fit by
#' iteratively reducing each EDF component while the penalized log-likelihood
#' improves by more than 2 units.
#'
#' @param y numeric response vector
#' @param x numeric predictor vector
#' @param startedf length-3 numeric vector; starting EDF values
#' @param boundedf length-3 numeric vector; minimum EDF values (search stops
#'   when an EDF component reaches its bound)
#' @param search.seq integer vector; order in which EDF components are searched
#'   (e.g. `c(1,3,2)` searches lambda, then sigma, then mu)
#' @param ... further arguments passed to [lmsqreg.fit()]
#' @return A list with components:
#'   \describe{
#'     \item{fit}{the best-fitting `lmsqreg.fit` object found}
#'     \item{pl}{its penalized log-likelihood}
#'     \item{edf}{the optimal EDF triple}
#'   }
#' @export
lmsqreg.search <- function(y, x, startedf = c(3, 5, 3),
    boundedf = c(1, 2, 1), search.seq = c(1, 3, 2), ...)
{
    oldfit <- lmsqreg.fit(y, x, edf = startedf, ...)
    oldpl  <- oldfit[[1]]$pl
    message(oldpl)
    newedf <- curedf <- startedf
    newpl  <- oldpl
    for (edfcomp in search.seq) {
        while (TRUE) {
            newedf[edfcomp] <- curedf[edfcomp] - 1
            if (newedf[edfcomp] < boundedf[edfcomp]) {
                newedf <- curedf
                break
            }
            message(paste(newedf, collapse = " "))
            tmpfit <- lmsqreg.fit(y, x, edf = newedf, ...)
            newpl  <- tmpfit[[1]]$pl
            message(newpl)
            if (newpl < (oldpl - 2))
                break
            else {
                oldfit <- tmpfit
                oldpl  <- newpl
                curedf <- newedf
            }
        }
    }
    list(fit = oldfit, pl = oldpl, edf = curedf)
}

#' Compute winsorization envelope bands
#'
#' Returns upper and lower boundary curves that enclose the non-outlier data
#' within each quantile-based section. Useful for overlaying winsorization
#' bands on a scatter plot.
#'
#' @param x numeric predictor vector
#' @param y numeric response vector
#' @param ncut integer; number of quantile-based sections
#' @param k integer; maximum number of outliers per section
#' @return A list with components `x`, `y` (upper band), `y2` (lower band),
#'   and `bad` (indices of outlying observations).
#' @export
locwin.envelope <- function(x, y, ncut = 5, k = 20)
{
    lamtab <- function(n, k, alpha = 0.05) {
        l <- 0:(k - 1)
        p <- 1 - ((alpha/2)/(n - l))
        tvals <- qt(df = n - l - 2, p = p)
        lams <- ((n - l - 1) * tvals)/sqrt((n - l - 2 + tvals^2) * (n - l))
        lams
    }
    skesd <- function(x, k) {
        bigres <- rep(NA, k)
        bigresind <- rep(NA, k)
        n <- length(x)
        curx <- x
        ind <- 1:n
        inds <- NULL
        for (i in 1:k) {
            last <- n - i + 1
            m <- mean(curx)
            ares <- abs(curx - m)
            oares <- order(ares)
            bigres[i] <- ares[oares[last]]/sqrt(var(curx))
            curx <- curx[-oares[last]]
            inds <- c(inds, ind[oares[last]])
            ind  <- ind[-oares[last]]
        }
        if (length(inds) == 0) {
            inds <- NA
            res  <- NA
        }
        list(res = bigres, ind = inds)
    }
    sgesdri <- function(x, k = ((length(x) %% 2) * floor(length(x)/2) +
        (1 - (length(x) %% 2)) * (length(x)/2 - 1)),
        alpha = 0.05)
    {
        E <- skesd(x, k)
        R <- E$res
        I <- E$ind
        n <- length(x)
        L <- lamtab(length(x), k, alpha = alpha)
        worst <- max((1:k)[R > L])
        if (is.na(worst) | is.null(worst))
            list(ind = NA, val = NA)
        else list(ind = I[1:worst], val = x[I[1:worst]])
    }
    cutfs <- (0:ncut)/ncut
    lims <- quantile(x, cutfs)
    lims[1] <- lims[1] - 0.1
    lims[length(lims)] <- lims[length(lims)] + 0.1
    xc <- cut(x, lims)
    iinds <- 1:length(x)
    iindsspl <- split(iinds, xc)
    xspl  <- split(x, xc)
    yspl  <- split(y, xc)
    yspl2 <- list()
    nout  <- list()
    bad   <- list()
    for (i in 1:length(yspl)) {
        curo <- sgesdri(yspl[[i]])$ind
        bad[[i]] <- iindsspl[[i]][curo]
        if (!is.na(curo[1])) {
            nout[[i]] <- length(curo)
            curcln    <- yspl[[i]][-curo]
        } else {
            nout[[i]] <- 0
            curcln    <- yspl[[i]]
        }
        lowcln   <- min(curcln)
        hicln    <- max(curcln)
        yspl[[i]]  <- rep(hicln,  length(yspl[[i]]))
        yspl2[[i]] <- rep(lowcln, length(yspl[[i]]))
    }
    message(paste(sum(unlist(nout)), "outliers"))
    retord <- order(as.integer(unlist(iindsspl)))
    list(
        x   = as.double(unlist(xspl))[retord],
        y   = as.double(unlist(yspl))[retord],
        y2  = as.double(unlist(yspl2))[retord],
        bad = unlist(bad)
    )
}

#' Validate an LMS fit on held-out data
#'
#' Tests goodness-of-fit of an `lmsqreg.fit` object on new data by comparing
#' nominal and actual quantile coverage and running KS, t-, and unit-variance
#' chi-squared tests within quantile-based x intervals.
#'
#' @param qreg an object of class `"lmsqreg.fit"`
#' @param y.val numeric vector of new response values
#' @param t.val numeric vector of new predictor values (same length as `y.val`)
#' @param rule integer; extrapolation rule passed to [stats::approx()]
#' @param xcuts numeric vector of quantile cut points for interval tests
#' @return Invisibly returns a list with quantile coverage (`pout`, `p.val`)
#'   and test p-value vectors `ps`, `tps`, `vps`.
#' @export
validate.report <- function(qreg, y.val, t.val,
    rule = 2, xcuts = quantile(t.val, c(0.2, 0.4, 0.6, 0.8)))
{
    qsys <- qreg$qsys
    pout <- rep(NA, length(qsys$pcts))
    for (k in 1:length(qsys$pcts)) {
        pk <- approx(qsys$targetx, qsys$outmat[k, ], t.val, rule = rule)
        pout[k] <- (sum(y.val < pk$y))/length(y.val)
    }
    z  <- zscores(y.val, t.val, qreg)
    TT <- t.val
    xcuts <- quantile(t.val, c(0.2, 0.4, 0.6, 0.8))
    xfac  <- cut(TT, round(c(min(TT) - 0.001, xcuts, max(TT) + 0.001), 3))
    zspl  <- split(z, xfac)
    ntests <- length(zspl)
    ps  <- rep(NA, ntests + 1)
    tps <- rep(NA, ntests + 1)
    vps <- rep(NA, ntests + 1)
    names(ps)  <- c(names(table(xfac)), "Overall")
    names(tps) <- names(ps)
    names(vps) <- names(ps)
    unit.var.test <- function(x, nullv = 1) {
        V <- var(x)
        n <- length(x)
        if (V <= nullv)
            ans <- (2 * pchisq(((n - 1) * V)/nullv, n - 1))
        else
            ans <- 2 * (1 - pchisq(((n - 1) * V)/nullv, n - 1))
        if (ans > 1) ans <- 1
        list(var = V, p.val = ans)
    }
    for (i in 1:ntests) {
        ps[i]  <- ks.test(zspl[[i]], "pnorm")$p.val
        tps[i] <- t.test(zspl[[i]])$p.val
        vps[i] <- unit.var.test(zspl[[i]])$p.val
    }
    ps[ntests + 1]  <- ks.test(z, "pnorm")$p.val
    tps[ntests + 1] <- t.test(z)$p.val
    vps[ntests + 1] <- unit.var.test(z)$p.val
    cat("Nominal quantile coverage\n")
    cat(round(qsys$pcts, 5))
    cat("\n")
    cat("Actual quantile coverage\n")
    cat(round(pout, 5))
    cat("\n")
    cat("\nKS tests: (intervals in", qreg[[1]]$xname, "//p-values)\n")
    print(round(ps, 3))
    cat("\nt tests: (intervals in", qreg[[1]]$xname, "//p-values)\n")
    print(round(tps, 3))
    cat("\nX2 tests (unit variance): (intervals in", qreg[[1]]$xname, "//p-values)\n")
    print(round(vps, 3))
    invisible(list(pout, p.val = qsys$pcts, ps = ps, tps = tps, vps = vps))
}

# Internal: extract the "version" attribute set on a function object.
Version <- function(x) attr(x, "version")
