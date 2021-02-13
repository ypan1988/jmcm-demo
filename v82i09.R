library("jmcm")
library("lattice")

##################################################
########## EXAMPLE: CATTLEDATA ####################
##################################################

data("cattle", package = "jmcm")
head(cattle)

xyplot(weight ~ day | group, group = id, data = cattle,
       xlab = "days", ylab = "weight", col = 1, type = "l")
cattleA <- subset(cattle, group == "A")
?jmcm

# mcd
# 1, 2, ..., 10, 10.5
system.time(fit1 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA, 
                         triple = c(8, 3, 4), cov.method = "mcd",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
# without output for the optimization
fit1 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
             triple = c(8, 3, 4), cov.method = "mcd")
# full likelihood
fit1a <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
              triple = c(8, 3, 4), cov.method = "mcd",
              control = jmcmControl(ignore.const.term = FALSE))
fit1
print(fit1) # output the same result as fit1
str(fit1)   # display the structure of fit1

getJMCM(fit1, "m")
getJMCM(fit1, "Y")
getJMCM(fit1, "X")
getJMCM(fit1, "Z")
getJMCM(fit1, "W")

# get the corresponding part for 15th subject
getJMCM(fit1, "m", 15)
getJMCM(fit1, "Y", 15)
getJMCM(fit1, "X", 15)
getJMCM(fit1, "Z", 15)
getJMCM(fit1, "W", 15)

getJMCM(fit1, "triple")

getJMCM(fit1, "theta")
getJMCM(fit1, "beta")
getJMCM(fit1, "lambda")
getJMCM(fit1, "gamma")

# get the estimated mean and covariance matrix for the 1st subject
getJMCM(fit1, "mu", 1)
getJMCM(fit1, "Sigma", 1)

# get the estimated mean and covariance matrix for the 30th subject
getJMCM(fit1, "mu", 30)
getJMCM(fit1, "Sigma", 30)

# meanplot(fit1)
regressogram(fit1, 1:11)

# acd
# 1, 2, ..., 10, 10.5
system.time(fit2 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
                         triple = c(8, 3, 4), cov.method = "acd",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
fit2 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
             triple = c(8, 3, 4), cov.method = "acd")
fit2
# meanplot(fit2)
regressogram(fit2, 1:11)

# hpc
# 1, 2, ..., 10, 10.5
system.time(fit3 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
                         triple = c(8, 2, 2), cov.method = "hpc",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
fit3 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
             triple = c(8, 2, 2), cov.method = "hpc")
fit3
# meanplot(fit3)
regressogram(fit3, 1:11)

##################################################
##########    TABLE 1   ##########################
##################################################

triples <- list(c(8, 3, 4),
                c(8, 2, 2),
                c(10, 10, 10),
                c(6, 1, 1),
                c(3, 3, 3),
                c(4, 4, 3),
                c(7, 2, 2),
                c(8, 7, 4),
                c(9, 1, 3),
                c(9, 4, 3),
                c(9, 8, 5))

table1 <- lapply(triples, function(triple) {
    results <- sapply(c("mcd", "acd", "hpc"), function(cov.method) {
        time <- system.time(fit <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA,
                                        triple = triple, cov.method = cov.method,
                                        control = jmcmControl(trace = TRUE, profile = TRUE)))
        list(fit = fit, time = time["elapsed"])
    }, simplify = FALSE)
    c(setNames(triple, c("p", "d", "q")),
      "No. of parsm" = sum(triple) + 3,
      unlist(lapply(results, function(x)
          c(l = x$fit@opt$loglik, BIC = x$fit@opt$BIC, time = x$time))))
})
do.call("rbind", table1)

##################################################
##########    EXAMPLE: CD4+   ####################
##################################################
data("aids", package = "jmcm")
head(aids, 20)

xyplot(sqrt(cd4) ~ time, data = aids,
  panel = function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.lines(x[aids$id == 10002], y[aids$id == 10002], col = 2, lwd = 2)
    panel.lines(x[aids$id == 10005], y[aids$id == 10005], col = 3, lwd = 2)
    panel.lines(x[aids$id == 10029], y[aids$id == 10029], col = 4, lwd = 2)
    panel.lines(x[aids$id == 10039], y[aids$id == 10039], col = 5, lwd = 2)
    panel.lines(x[aids$id == 10048], y[aids$id == 10048], col = 6, lwd = 2)
    panel.lines(x[aids$id == 10052], y[aids$id == 10052], col = 7, lwd = 2)
  },
  xlab = "Time", ylab = "CD4 cell numbers", col = 1)

# mcd
system.time(fit4 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
                         triple = c(8, 1, 3), cov.method = "mcd",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
# without output for the optimization
fit4 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
             triple = c(8, 1, 3), cov.method = "mcd")
fit4

fit4a <- jmcm(I(sqrt(cd4)) | id | time ~ age | age + packs, data = aids,
              triple = c(8, 1, 3), cov.method = "mcd")

# check the difference between the covariate matrix for the 1st subject 
getJMCM(fit4, "X", 1)
getJMCM(fit4a, "X", 1)
getJMCM(fit4, "Z", 1)
getJMCM(fit4a, "Z", 1)

# regressogram(fit4)
meanplot(fit4)
bootcurve(fit4, nboot = 1000) 

# acd
system.time(fit5 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
                         triple = c(8, 1, 3), cov.method = "acd",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
fit5 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
             triple = c(8, 1, 3), cov.method = "acd")
fit5

# regressogram(fit5)
meanplot(fit5)
bootcurve(fit5, nboot = 1000) # it may take quite a long time

# hpc
system.time(fit6 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
                         triple = c(8, 1, 1), cov.method = "hpc",
                         control = jmcmControl(trace = TRUE, profile = TRUE)))
fit6 <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
             triple = c(8, 1, 1), cov.method = "hpc")
fit6

# regressogram(fit6)
meanplot(fit6)
bootcurve(fit6, nboot = 1000) # it may take hours

##################################################
##########    TABLE 2   ##########################
##################################################

triples <- list(c(8, 1, 1),
                c(8, 1, 3),
                c(6, 1, 1),
                c(3, 3, 3),
                c(4, 4, 3),
                c(8, 3, 3),
                c(8, 7, 4),
                c(9, 1, 3),
                c(9, 4, 3),
                c(9, 8, 5))

table2 <- lapply(triples, function(triple) {
    results <- sapply(c("mcd", "acd", "hpc"), function(cov.method) {
        time <- system.time(fit <- jmcm(I(sqrt(cd4)) | id | time ~ 1 | 1, data = aids,
                                        triple = triple, cov.method = cov.method,
                                        control = jmcmControl(trace = TRUE, profile = TRUE)))
        list(fit = fit, time = time["elapsed"])
    }, simplify = FALSE)
    c(setNames(triple, c("p", "d", "q")),
      "No. of parsm" = sum(triple) + 3,
      unlist(lapply(results, function(x)
          c(l = x$fit@opt$loglik, BIC = x$fit@opt$BIC, time = x$time))))
})
do.call("rbind", table2)


