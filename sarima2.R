# sarima2
# Description: Slight rewrite of astsa::sarima() to modularize the code. The 
# plot and Ljung-Box statistics are split from the main code, allowing use 
# beyond the sarima() function. Additional arguments added to control the fit 
# report, control the plot, and offset the trend. Some variables were renamed
# in order to clarify meaning and compare with Python's statsmodels ARIMA 
# function. The plot function, plot.Arima, is also extended to include a 
# residual histogram and the PACF of the residuals. It can be used on some non-
# ARIMA models, provided the sigma() function works.

# Contents:
# - read.ts - Unnecessary function to read in time series quietly.
# - get_lags - Helper function to fetch lags, not intended for general use
# - sarima_title - Helper function to assemble ARIMA title, not intended for general use
# - ljung_box - Ljung-Box test separated from sarima() so it can be used elsewhere
# - plot.Arima - Diagnostic plots separated from sarima(), which can be used on
#   Arima model fits with plot() or non-Arima models with plot.Arima(). 
#   Automatic identification of model parameters can be overridden with the
#   parameters= argument.
# - sarima2 - Modularized sarima(), backward compatible.

# Function to *quietly* read in time series data
# Note that all additional arguments are sent to ts() and not scan()
read.ts <- function(input_string, ...) {
  output_object <- ts(scan(input_string, quiet=T), ...)
  output_object
}

# Helper function to calculate lags. Not intended for general use.
get_lags <- function(fitit, lags=NULL, nfixed=0) {
  
  num <- sum(!is.na(residuals(fitit)))
  
  # Set lags if not specified
  if ("Arima" %in% class(fitit)) {
    S <- fitit$arma[5]
    # ppq <- p + q + P + Q - nfixed
    ppq <- fitit$arma[1] + fitit$arma[2] + fitit$arma[3] + fitit$arma[4] - nfixed
    if (is.null(lags)){
      nlag <- ifelse(S < 7, 20, 3 * S)
      nlag <- ifelse(nlag < ppq + 8, ppq + 8, nlag)
    } else {
      nlag = lags
    }
  } else {
    S <- NULL
    if (is.null(lags)){
      alag <- 10 + sqrt(num)
      nlag <- alag
    } else {
      alag <- lags
      nlag <- lags
    }
  }
  
  # Autocorrelation lags
  # Assess lags and override with user-specified setting
  if (is.null(lags)){
    alag <- max(10 + sqrt(num), 3 * S)
  } else {
    alag <- lags
  }
  
  c(nlag, alag)
}

# Assemble the ARIMA model title: "(p, d, q)(P, D, Q) [S]"
sarima_title <- function(fitit) {
  p <- fitit$arma[1]
  q <- fitit$arma[2]
  P <- fitit$arma[3]
  Q <- fitit$arma[4]
  S <- fitit$arma[5]
  d <- fitit$arma[6]
  D <- fitit$arma[7]
  if (S <= 0) {
    mtitle <- paste0(
      "Model: (", p, ",", d, ",", q, ")"
    )
  } else {
    mtitle <- paste0(
      "Model: (", p, ",", d, ",", q, ") ", 
      "(", P, ",", D, ",", Q, ") [", S, "]"
    )
  }
}

# Ljung-Box Statistics, provided as a data frame
ljung_box <- function(fitit, lags=NULL, nfixed=0) {
  
  if ('Arima' %in% class(fitit)){
    ppq <- fitit$arma[1] + fitit$arma[2] + fitit$arma[3] + fitit$arma[4] - nfixed
  } else {
    parameters <- length(fitit$coefficients[names(fitit$coefficients) != '(Intercept)'])
    ppq <- parameters - nfixed
  }
  
  nlag <- get_lags(fitit, lags=lags, nfixed=nfixed)[1]
  qval <- numeric(nlag)
  pval <- numeric(nlag)
  
  for (i in (ppq + 1):nlag) {
    qval[i] <- stats::Box.test(residuals(fitit), i, type = "Ljung-Box")$statistic
    pval[i] <- stats::pchisq(qval[i], i - ppq, lower.tail = FALSE)
  }
  
  data.frame('qval'=qval, 'pval'=pval)
  
}

# Diagnostic plots for assessing Arima objects.
# Can be used on non-arima models by specifying plot.Arima(fitted_object)
plot.Arima <- function(fitit, lags=NULL, nfixed=0, S=NULL, parameters=NULL, main=NULL, extended=TRUE, ...) {

  # Residuals and number of observations
  rs <- residuals(fitit)
  num <- sum(!is.na(rs))
  lag_list <- get_lags(fitit, lags=lags, nfixed=nfixed)
  nlag <- lag_list[1]
  alag <- lag_list[2]

  # for ARIMAs, extract specifications
  if ("Arima" %in% class(fitit)) {

    # standardized residuals
    stdres <- rs / sqrt(fitit$sigma2)

    # $arma specifies AR, MA, SAR, SMA, S, diff, Sdiff
    p <- fitit$arma[1]
    q <- fitit$arma[2]
    P <- fitit$arma[3]
    Q <- fitit$arma[4]
    S <- fitit$arma[5]
    d <- fitit$arma[6]
    D <- fitit$arma[7]
    ppq <- p + q + P + Q - nfixed
    num <- fitit$nobs
    
    if (is.null(main)){
      # Create a main title with the model specifications
      main <- sarima_title(fitit)
    }
        
  } else {
    
    parameters <- ifelse(
      is.null(parameters), 
      length(fitit$coefficients[names(fitit$coefficients) != '(Intercept)']), 
      parameters
    )
    ppq <- parameters - nfixed
    num <- length(residuals(fitit))
    stdres <- residuals(fitit) / fitit$sigma2

  }

  layout = graphics::layout
  par = graphics::par
  plot = graphics::plot
  acf = stats::acf
  polygon = graphics::polygon
  abline = graphics::abline
  lines = graphics::lines
  title = graphics::title
  frequency = stats::frequency
  coef = stats::coef
  dnorm = stats::dnorm
  ppoints = stats::ppoints
  qnorm = stats::qnorm
  time = stats::time
  na.pass = stats::na.pass
  
  old.par <- par(no.readonly = TRUE)
  if (extended) {
    layout(matrix(c(1, 2, 4, 6, 1, 3, 5, 6), ncol = 2))
  } else {
    layout(matrix(c(1, 2, 4, 1, 3, 4), ncol = 2))
  }

  # 1. Residuals time series plot
  tsplot(stdres, main = "Standardized Residuals", ylab = "", ...)
  title(main, adj = 0)
  title('plot.Arima', adj=1)

  # 2. ACF plot
  acf1(rs, alag, main = "ACF of Residuals", ...)
  
  # 3. Q-Q plot
  u = qqnorm(stdres, plot.it = FALSE)
  lwr = min(-4, min(stdres))
  upr = max(4, max(stdres))
  tsplot(
    u$x, u$y, type = "p", ylim = c(lwr, upr), 
    ylab = "Sample Quantiles", 
    xlab = "Theoretical Quantiles", 
    main = "Normal Q-Q Plot of Std Residuals", 
    ...
  )
  # Q-Q plot diagonal line
  sR <- !is.na(stdres)
  ord <- order(stdres[sR])
  ord.stdres <- stdres[sR][ord]
  PP <- stats::ppoints(num)
  z <- stats::qnorm(PP)
  y <- stats::quantile(
    ord.stdres, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE
  )
  x <- stats::qnorm(c(0.25, 0.75))
  b <- diff(y)/diff(x)
  a <- y[1L] - b * x[1L]
  abline(a, b, col = 4)
  # Q-Q plot confidence band
  SE <- (b/dnorm(z)) * sqrt(PP * (1 - PP)/num)
  qqfit <- a + b * z
  U <- qqfit + 3.9 * SE
  L <- qqfit - 3.9 * SE
  z[1] = z[1] - 0.1
  z[length(z)] = z[length(z)] + 0.1
  xx <- c(z, rev(z))
  yy <- c(L, rev(U))
  polygon(xx, yy, border = NA, col = gray(0.5, alpha = 0.2))
  
  # 4. PACF plot
  if (extended) {
    acf1(rs, alag, main = "PACF of Residuals", pacf = TRUE, ...) 
  }
  
  # 5. Histogram
  if (extended) {
    
    if (length(stdres) > 5000) {
      spresults <- ''
    } else {
      sp1 <- shapiro.test(stdres)$`p.value`
      if (round(sp1, 4) == 0){
        spresults <- paste0(' (S-W p<0.001)')
      } else {
        spresults <- paste0(' (S-W p=', round(sp1, 4), ')')
      }
    }
    
    d1 <- density(stdres)
    hist(
      stdres, freq=F, 
      main = paste0('Histogram of Std Residuals', spresults), 
      xlab = 'Std. Residuals', ylim=c(0, max(d1$y))
    )
    lines(d1)
  }
  
  # 6. Ljung-Box p-value plot
  # This is plot 4 if extended=FALSE
  pval <- ljung_box(fitit, nfixed=nfixed)$pval
  tsplot(
    (ppq + 1):nlag, pval[(ppq + 1):nlag], type = "p",
    xlab = "LAG (H)", ylab = "p value", ylim = c(-0.14, 1), 
    main = "p values for Ljung-Box statistic", 
    ...
  )
  abline(h = 0.05, lty = 2, col = 4)
  
  on.exit(par(old.par))
  
}

sarima2 <- function(
    xdata, p, d, q, P = 0, D = 0, Q = 0, S = -1, 
    no.constant = FALSE, xreg = NULL, fixed = NULL, 
    Model = TRUE, plot = TRUE, report = FALSE, details = TRUE, extended = TRUE,
    trend_offset = 1,
    tol = sqrt(.Machine$double.eps), ...
) {
  
  arima = stats::arima
  trans = ifelse(is.null(fixed), TRUE, FALSE)
  n = length(xdata)

  plot = ifelse(details, plot, FALSE)
  report = ifelse(details, report, FALSE)
  trc = ifelse(report, 1, 0)

  # With no external regressors
  if (is.null(xreg)) {
    # trend is used as a linear time factor
    # trend used to be called constant
    trend = (trend_offset):(n + trend_offset - 1)
    # constant is used to estimate the mean
    # constant use to be called xmean
    constant = rep(1, n)
    if (no.constant == TRUE) {
      constant = NULL
    }
    
    # No differences
    # Uses a series of 1s as an external regressor to estimate mean
    if (d == 0 & D == 0) {
      fitit = arima(
        xdata, 
        order = c(p, d, q), 
        seasonal = list(order = c(P, D, Q), period = S), 
        xreg = constant, # estimate mean
        include.mean = FALSE, # set to false to control with no.constant
        # fixed and transform.pars are paired
        # if fixed is set, transform.pars is FALSE
        fixed = fixed, # fixed values of ARIMA coefficients
        transform.pars = trans, # transform AR parameters for stationarity
        optim.control = list(trace = trc, REPORT = 1, reltol = tol)
      )
    }
    # Either a nonseasonal difference or a seasonal difference, not both
    # Uses time indicators (e.g, 1, 2, 3, 4...) for intercept estimate
    else if (xor(d == 1, D == 1) & no.constant == FALSE) {
      fitit = arima(
        xdata, 
        order = c(p, d, q), 
        seasonal = list(order = c(P, D, Q), period = S), 
        xreg = trend, # estimate time factor
        #include.mean = FALSE, # this is ignored with differencing
        # fixed and transform.pars are paired
        # if fixed is set, transform.pars is FALSE
        fixed = fixed, # fixed values of ARIMA coefficients
        transform.pars = trans, # transform AR parameters for stationarity
        optim.control = list(trace = trc, REPORT = 1, reltol = tol)
      )
    }
    # Both differences
    # Despite the code below, does not generate mean or constant regardless of user settings
    else {
      fitit = arima(
        xdata, 
        order = c(p, d, q), 
        seasonal = list(order = c(P, D, Q), period = S), 
        include.mean = !no.constant, # this actually does nothing regardless of setting
        # fixed and transform.pars are paired
        # if fixed is set, transform.pars is FALSE
        fixed = fixed, # fixed values of ARIMA coefficients
        transform.pars = trans, # transform AR parameters for stationarity
        optim.control = list(trace = trc, REPORT = 1, reltol = tol)
      )
    }
  }
  # With external regressors
  # Does not estimate a mean, nor does it use a time indicator variable
  else {
    fitit = stats::arima(
      xdata, 
      order = c(p, d, q), 
      seasonal = list(order = c(P, D, Q), period = S), 
      xreg = xreg, 
      # fixed and transform.pars are paired
      # if fixed is set, transform.pars is FALSE
      fixed = fixed, # fixed values of ARIMA coefficients
      transform.pars = trans, # transform AR parameters for stationarity
      optim.control = list(trace = trc, REPORT = 1, reltol = tol)
    )
  }

  # Create a main title with the model specifications
  mtitle <- sarima_title(fitit)
  
  # Diagnostic plot
  if (plot) {
    # Generate the diagnostic plot
    nfixed = sum(!is.na(fixed))
    if (Model) {
      plot.Arima(fitit, nfixed=nfixed, main=mtitle, extended=extended, ...)
    } else {
      plot.Arima(fitit, nfixed=nfixed, extended=extended, ...)
    }
    
  }
  
  # t-test results
  if (is.null(fixed)) {
    coefs = fitit$coef
  }
  else {
    coefs = fitit$coef[is.na(fixed)]
  }
  k = length(coefs)
  n = fitit$nobs
  dfree = n - k
  t.value = coefs / sqrt(diag(fitit$var.coef))
  p.two = stats::pf(t.value^2, df1 = 1, df2 = dfree, lower.tail = FALSE)
  ttable = cbind(
    Estimate = coefs, SE = sqrt(diag(fitit$var.coef)), t.value, p.value = p.two
  )
  ttable = round(ttable, 4)
  
  # IC values
  BIC = stats::BIC(fitit)/n
  AIC = stats::AIC(fitit)/n
  AICc = (n * AIC + ((2 * k^2 + 2 * k)/(n - k - 1)))/n
  
  # Final output
  list(
    fit = fitit, degrees_of_freedom = dfree, ttable = ttable, 
    AIC = AIC, AICc = AICc, BIC = BIC,
    title = mtitle
  )

}

sarima2.for <- function(xdata, n.ahead, p=NA, d=NA, q=NA, P=0, D=0, Q=0, S=-1,
  modelObject=NA, tol=sqrt(.Machine$double.eps), 
  no.constant=FALSE, plot=TRUE, plot.all=FALSE,
  xreg=NULL, newxreg=NULL, fixed=NULL, ...
) { 
  
  #xdata <- bh.train$Hires2; modelObject <- mod4ts
  
  if (is.na(p)) {
    p <- modelObject$fit$arma[1]
    q <- modelObject$fit$arma[2]
    P <- modelObject$fit$arma[3]
    Q <- modelObject$fit$arma[4]
    S <- modelObject$fit$arma[5]
    d <- modelObject$fit$arma[6]
    D <- modelObject$fit$arma[7]
  }
  
  sarima.for(
    xdata=xdata, 
    n.ahead=n.ahead, p=p, d=d, q=q, P=P, D=D, Q=Q, S=S, 
    tol=tol, 
    no.constant=no.constant, plot=plot, plot.all=plot.all, 
    xreg=xreg, newxreg=newxreg, fixed=fixed, 
    ...
  )

}
