### Biostatistics

gamma(100)
?formula


## simple linear regression
fit <- lm(tm ~ gc, data = md)
predict(fit, data.frame(gc=42))

function(x, pred, linecol="darkred", fillcol="gold", alpha=0.5) {
  # fill the area between lower and upper
  # the filling colour can be made transparent
  # but alpha must be scaled from 0..1 to 0..255 as col2rgb provides values in that range
  # and what it returns must be transposed, too
  fcol <- rgb(t(col2rgb(fillcol)), max=255, alpha=255*alpha)
  polygon(c(x,rev(x)), c(pred[,2],rev(pred[,3])), border=NA, col=fcol)
  
  # these are the predicted values
  lines(x, pred[,1], lty="solid", col=linecol, lwd=2)
  # lower and upper range, usually confidence intervals
  lines(x, pred[,2], lty="dashed", col=linecol)
  lines(x, pred[,3], lty="dashed", col=linecol)
}

## Linear model with two 
lm(tm ~ gc + lna, data = md2)

## linear model with interactions
lm(tm ~ gc * lna, data = md2)

## oNe can check the likelihood for these models using AIC or BIC

## Fit non-linear functions when the parameters are linear (a, b --- coefficients)
lm(y ~ I(sin(pi*x/12)))  ## fitting sine function in R


## Polynomials
x <- 0:96
p <- sapply(1:4, function(n){x^n})
cor(p)  ## you can see very high correlations

## orthogonal polynomials
op  <- poly(x, degree = 4)
cor(op)


##we do orthogonal polynomial fit
?poly

##__________________________________________
# ANOVA

#Group mean and sample mean
## To test for the homoscedasticity test
bartlett.test(values ~ ind, sc)
fligner.test(values ~ ind, sc)

dd_aov <- aov(values ~ ind, sc)
TukeyHSD(dd_aov)
plot(TukeyHSD(dd_aov))

### Weighted linear regression


### Special case t-test
set.seed(137)
g1 <- rnorm(20, mean = 1.7, sd = 0.4)
g2 <- rnorm(20, mean = 2.1, sd = 0.4)
wf <- data.frame(grp1=g1, grp2 = g2)
tf <- stack(wf)

t.test(g1, g2, var.equal = T)
summary(aov(values ~ ind, tf))

## what if they have different sd
g2 <- rnorm(20, mean = 2.1, sd = 0.6)


cols <- c("green", "orange", "blue")
### Special case t-test
set.seed(137)
g1 <- rnorm(20, mean = 1.7, sd = 0.4)
g2 <- rnorm(20, mean = 2.1, sd = 0.4)
wf <- data.frame(grp1=g1, grp2 = g2)
tf <- stack(wf)

t.test(g1, g2, var.equal = T)
summary(aov(values ~ ind, tf))

## what if they have different sd
g2 <- rnorm(20, mean = 2.1, sd = 0.6)

###Try me! unequal variances
df <- read.csv("heterosced.csv")
boxplot(df, col = cols)

bartlett.test(df)
fligner.test(df)

### Welch's one-way test, similar to ANOVA with different variances
sf <- stack(df)
print(oneway.test(values ~ ind, sf))

###Two way ANOVA
d2 <- read.csv("a2.csv")
ms <- with(d2, tapply(growth, list(genotype, feed), mean))
barplot(ms, beside = T, col = cols)

model <- aov(growth ~ genotype * feed, d2)
summary.aov(model)
summary.lm(model)

noint.model <- aov(terms(growth ~ feed + genotype, keep.order = T), d2)
summary(aov(growth ~ genotype + feed, d2))
summary.aov(noint.model)

hsd <- TukeyHSD(noint.model, "genotype")
print(hsd)
plot(hsd)

######Generalised linear models
df <- read.csv("deadmice.csv")
head(df)
d <- t(as.matrix(df[,2:3]))
barplot(d, beside = T, col = cols[2:3])
plot(df$dose, df$prop, col = cols[3], pch = 19)

# fitting a binomial 

######Generalised linear models
df <- read.csv("deadmice.csv")
head(df)
d <- t(as.matrix(df[,2:3]))
barplot(d, beside = T, col = cols[2:3])
plot(df$dose, df$prop, col = cols[3], pch = 19)

# fitting a binomial GLM 
x <- df$dose
y <- cbind(df$killed, df$treated - df$killed)
gfit <- glm(y ~ x, family = binomial)
summary(gfit)


df <- read.csv("deadmice.csv")
head(df)
d <- t(as.matrix(df[,2:3]))
barplot(d, beside = T, col = cols[2:3])
plot(df$dose, df$prop, col = cols[3], pch = 19)

# fitting a binomial GLM 
x <- df$dose
y <- cbind(df$killed, df$treated - df$killed)
gfit <- glm(y ~ x, family = binomial)
summary(gfit)

# dont understand the difference
plot(x, y[,1])
abline(lm(y[,1] ~ x))

##Counting reads
df <- read.csv("poissoncount.csv")
plot(df, pch = 19, col = "blue")

# Quick fix -1, take log of the values
log.model <- lm(log(y) ~ x, data = df)
summary(log.model)
AIC(log.model)

# quicj fix-2, weighted based on the weights 1/y
wgt.model <- lm(y ~ x, data = df, weights = 1/y)
summary(wgt.model)
AIC(wgt.model)

# poisson glm
glm.model <- glm(y ~ x, data = df, family = poisson)
glm.model <- glm(y ~ x, data = df, family = gaussian)  ## taking an error model as gaussian, 
#you can see very high residual deviance
summary(glm.model)

### Counting real-world reads!
df <- read.csv("overdispcount.csv")
plot(df, pch = 19, col = "blue")

glm.model <- glm(y ~ x, data = df, family = poisson)
glm.model <- glm(y ~ x, data = df, family = quasipoisson)
summary(glm.model)
## negative binomial glm, better distribution for over dispersion poisson
library("MASS")
nb.model <- glm.nb(y ~ x, data = df)
summary(nb.model)

#_________________________________
#### ANCOVA
# a hybrid between the two variables
diff <- read.csv("ko_diff.csv")

fit.wt <- lm(y ~ x, data = diff, subset = which(ind == "wt"))
summary(fit.wt)

fit.ko <- lm(y ~ x, data = diff, subset = which(ind == "ko"))
summary(fit.ko)

fit <- lm(y ~ x * ind, data = diff)
summary(fit)
anova(fit)

source("ancovamodels.R")
print.ancova("ko_slope.csv")



