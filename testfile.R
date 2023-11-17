library(dmhBGM)
data("Wenchuan")
# nrow(na.omit(Wenchuan))

t0 <- proc.time()
a <- dmhbgm(x = Wenchuan, dmhsamples = 20, iter = 10, burnin = 1, parallel = FALSE)
t1 <- proc.time()
(t1 - t0) / 60

RcppParallel::defaultNumThreads()
RcppParallel::setThreadOptions(8)

t0 <- proc.time()
b <- dmhbgm(x = Wenchuan, dmhsamples = 20, iter = 10, burnin = 1, parallel = TRUE)
t1 <- proc.time()
(t1 - t0) / 60
