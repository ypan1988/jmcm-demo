library("jmcm")

cat("\n\n   Jmcm Benchmark 1.0\n")
cat("   ==================\n")
cat("\n\n")

cat("   I. Cattle Data\n")
cat("   ---------------------\n")
if (R.Version()$os == "Win32" || R.Version()$os == "mingw32") flush.console()

data("cattle", package = "jmcm")
cattleA <- subset(cattle, group == "A")
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
                                        control = jmcmControl(trace = FALSE, profile = TRUE)))
        list(fit = fit, time = time["elapsed"])
    }, simplify = FALSE)
    c(setNames(triple, c("p", "d", "q")),
      "No. of parsm" = sum(triple) + 3,
      unlist(lapply(results, function(x)
          c(l = x$fit@opt$loglik, BIC = x$fit@opt$BIC, time = x$time))))
})
print(do.call("rbind", table1))

cat("\n")

cat("   II. CD4+ Data\n")
cat("   ---------------------\n")
if (R.Version()$os == "Win32" || R.Version()$os == "mingw32") flush.console()

data("aids", package = "jmcm")
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
                                        control = jmcmControl(trace = FALSE, profile = TRUE)))
        list(fit = fit, time = time["elapsed"])
    }, simplify = FALSE)
    c(setNames(triple, c("p", "d", "q")),
      "No. of parsm" = sum(triple) + 3,
      unlist(lapply(results, function(x)
          c(l = x$fit@opt$loglik, BIC = x$fit@opt$BIC, time = x$time))))
})
print(do.call("rbind", table2))

cat("\n                      --- End of test ---\n\n")  
