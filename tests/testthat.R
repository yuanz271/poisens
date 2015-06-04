library(testthat)
library(poisens)
library(microbenchmark)
library(MASS)
library(glmnet)
library(plyr)
library(foreach)

# test_check("poisens2")

p <- read.csv("http://www.ats.ucla.edu/stat/data/poisson_sim.csv")
m1 <- glm(num_awards ~ prog + math, family="poisson", data=p)
coef(m1)
predict(m1, type = "response")
x <- model.matrix(m1)
y <- m1$y
m2 <- glm(num_awards ~ prog + math, family=negative.binomial(1), data=p)
coef(m2)

(f1 <- pois_reg(x, y))
mean((predict(f1, x) - y)^2)
negbin_reg(x, y, 1)

e1 <- pois_enet(x, y, lambda = pois_lambda(x, y))
predict(e1, x)

ens1 <- pois_reg_ens(x, y)
mean((predict(ens1, x) - y)^2)
ens2 <- negbin_reg_ens(x, y, 1)
mean((predict(ens2, x) - y)^2)


m3 <- glmnet(x[,-1], y, family="poisson", nlambda = 10, alpha = 1, standardize = TRUE)
coef(m3)

pois_enet(x, y, alpha = 1, lambda = pois_lambda(x, y, 1, 10), intercept = TRUE, 5000, 1e-6)
negbin_enet(x, y, theta = 1, alpha = 1, lambda = negbin_lambda(x, y, theta = 1, 1, 10), intercept = T, 5000, 1e-6)

microbenchmark(predict(f1, x))
microbenchmark(exp(x %*% f1$coef))

ens3 <- pois_enet_ens(x, y, 1, pois_lambda(x, y, 1, 10))
predict(ens3, x)

library(MASS)
quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
x <- model.matrix(quine.nb1$terms, quine.nb1$model)
y <- quine.nb1$y
fit2 <- negbin_reg2(x,y)
fit2$phi_se

t1 <- quine.nb1$theta^2  / quine.nb1$SE.theta^2
t2 <- (fit2$phi / fit2$phi_se)^2
1 - pchisq(t2, 1)

counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
glm.D93 <- glm.nb(counts ~ outcome + treatment)
x <- model.matrix(glm.D93$terms, glm.D93$model)
y <- glm.D93$y
fit2 <- negbin_reg2(x,y)
fit3 <- negbin_reg_ens1(x ,y)
fit4 <- negbin_reg_ens2(x ,y)
1 - pchisq(as.numeric((fit2$phi/fit2$phi_se)^2), 1)

