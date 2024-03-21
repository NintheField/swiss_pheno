

gaussian_wrap <- function(dat,sigma){

  naind <- c(which(is.na(dat)),length(dat) +1)
  gausm_all <- rep(NA, length(dat))
  start <- 1

  for(ii in 1:length(naind)){
    sind <- start:(naind[ii]-1)
    if(!length(sind) < 5){gausm_all[sind] <- mmand::gaussianSmooth(dat[sind], sigma = sigma)}
    start <- naind[ii] + 1
  }

  return(gausm_all)
}
