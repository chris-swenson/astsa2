This is a rewrite of an older version of sarima() from the astsa package. I have modularized the code, splitting the plotting from fitting, as well as some helper functions that overlap between the two.

The latest version of sarima() from astsa has added the ability to run the diagnostic plot without fitting a model, which overlaps with the functionality I provide here; however, the plot.Arima function extracted still functions well, including on non-Arima objects, and it provides a default plot() method for Arima objects.

# Usage:  

## sarima2

sarima2(  
&nbsp;&nbsp;&nbsp;&nbsp;xdata, p, d, q, P = 0, D = 0, Q = 0, S = -1,  
&nbsp;&nbsp;&nbsp;&nbsp;no.constant = FALSE, xreg = NULL, fixed = NULL,  
&nbsp;&nbsp;&nbsp;&nbsp;Model = TRUE, plot = TRUE, report = FALSE, details = TRUE, extended = TRUE,  
&nbsp;&nbsp;&nbsp;&nbsp;trend_offset = 1,  
&nbsp;&nbsp;&nbsp;&nbsp;tol = sqrt(.Machine$double.eps), ...  
)
- Just like sarima(), with these differences in arguments: 
- details=: Functions the same as sarima(), producing the plots after fitting.
- plot=: Boolean indicator to produce the plots, intended to be used independently of details=.
- report=: Boolean indicator to output the trace result from the arima() call.
- extended=: Boolean indicator to extend the plots, adding a residual histogram with KDE curve and the PACF of the residuals.
- trend_offset=: The initial value to use in the trend factor (called "constant" in the original code). Useful to set it at a different value if you're starting from some other point than 1.

## sarima2.for

sarima2.for(  
&nbsp;&nbsp;&nbsp;&nbsp;xdata, n.ahead, p=NA, d=NA, q=NA, P=0, D=0, Q=0, S=-1,  
&nbsp;&nbsp;&nbsp;&nbsp;modelObject=NA, tol=sqrt(.Machine$double.eps),  
&nbsp;&nbsp;&nbsp;&nbsp;no.constant=FALSE, plot=TRUE, plot.all=FALSE,  
&nbsp;&nbsp;&nbsp;&nbsp;xreg=NULL, newxreg=NULL, fixed=NULL, ...  
)
- Used just like sarmia.for() but:
- modelObject=: A fitted model object, used to extract all the parameters so you don't have to repeat them.
- e.g.:  
`arima_model_1 <- sarima2(input_data, ...)`  
`sarima2.for(input_data, arima_model_1$fit)`

## plot.Arima

plot.Arima(  
&nbsp;&nbsp;&nbsp;&nbsp;fitit, lags=NULL, nfixed=0, S=NULL, parameters=NULL, main=NULL, extended=TRUE, ...  
)
- fitit: If the input is a fitted ARIMA model, it collects all the information from the model object. Otherwise, it uses sigma().
- lags=: Override the automatic calculation of how many lags to display with a set number.
- nfixed=: The number of fixed model parameters.
- S=: The seasonal frequency if the input is not an ARIMA model.
- parameters=: The number of parameters in the model if the automatic calculation isn't correct.
  - The automatic calculation is the number of coefficients excluding the intercept
- main=: The title of the plot, like usual.
  - The automatic title is based on the Arima object parameters, varying on whether S= is populated
  - Use this if your model is not an Arima and you want a fitting title
- extended=: Boolean indicator to extend the plots, adding a residual histogram with KDE curve and the PACF of the residuals.
