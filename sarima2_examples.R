library(astsa)
library(nlme)
source('./sarima2/sarima2.R')

# Example 1: ARIMA
birth <- astsa::birth
model1 <- sarima2(birth, 1, 1, 1, 2, 1, 1, 12, plot=F)
plot(model1$fit)
plot(model1$fit, extended=F)
sarima(birth, 1, 1, 1, 2, 1, 1, 12) # For comparison

# Example 2: Linear
diff1and12 <- diff(diff(birth, 1), 12)
birth_lag1 <- lag(diff1and12, -1)
birth_lag12 <- lag(diff1and12, -12)
newbirth <- ts.intersect(diff1and12, birth_lag1, birth_lag12)
model2 <- lm(newbirth[,1] ~ newbirth[,2] + newbirth[,3])
plot.Arima(model2, S=12, main='Linear Pseudo-ARIMA') # specify S since it's not a part of the LM object

# Example 3: Polynomial
model3 <- lm(birth ~ poly(time(birth), 3))
pred3 <- predict(model3)
plot(birth, type='l')
lines(x=time(birth), y=pred3, type='l', col='red')
plot.Arima(model3, S=12, main='Polynomial Order 3') # specify S since it's not a part of the LM object

# Example 4: GLS with correlated errors
diff1and12 <- diff(diff(birth, 1), 12)
model4 <- gls(diff1and12 ~ time(diff1and12), correlation = corARMA(form = ~ 1, p=1, q=1))
# The following specifications are just a test, they are not correct
plot.Arima(model4, parameters=5, nfixed=2)
