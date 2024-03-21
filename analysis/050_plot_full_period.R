#### plot the recons based on model an periods

rm(list=ls())

setwd("/scratch3/nimfeld/wear/swiss_pheno")

#### packages
library(ggplot2)
library(dplyr)
#library(hrbrthemes)
library(viridis)
library(reshape2)

classes <- "1234_longseries_filled"
burnin <- 4000
iter <- 36000
date <- "2024-03-08"
pv <- "mprua65d"

#### directories
base_dir <- file.path(
  "data/calibrations",
  paste0("bt_run_",date,"_class",classes,"_",burnin,"burnin_",iter,"iter")
)

recon_dir <- file.path(base_dir,pv)

files <- list.files(file.path(base_dir,pv), pattern = "model_comparison", full.names = T)
mod_nams <- paste0("P",stringr::str_sub(files,-12,-5))
mod_outs <- lapply(files, readRDS)
names(mod_outs) <- mod_nams

### temperature data
stn_files <- list.files(file.path(base_dir,pv), pattern = "stn_list_train", full.names = T)
stn_lists <- lapply(stn_files, readRDS)
names(stn_lists) <- mod_nams

doymat <- readRDS(paste0(recon_dir,"/liestal_estimate_",pv,"_periods_model_",classes,".rds"))

mods <-  c("LIN","TT", "TTs", "PTT", "PTTs","M1","M1s","AT")
periods <- c("1961-1980","1971-1990","1981-2000","1991-2010","2000-2020")

cols <- c("#FF6347", "#FFD700", "#20B2AA", "#9932CC", "#00CED1", "#FFA07A")
lwd <- 1

#### all models together
png(paste0(recon_dir,"/figures/mcomp_allmods_20yrperiods_",Sys.Date(),".png"), width = 1500, res = 300, height = 2000, pointsize = 7)
par(mfrow = c(length(mods),1), mar = c(0.5,2,0.5,1), oma= c(3,1,1,1))

for (mm in 1:length(mods)){
  plot(1763:2020,doymat[,1,1], type = "l", ylim = c(90, 135), lwd = lwd, col = "white", xaxt = ifelse(mm == 7, "s","n"))
  for(i in 1:length(periods)) points(1763:2020,doymat[,i,mm], type = "l", lwd = lwd, col = cols[i])
  mtext(side = 3, mods[mm], line = -1.2)
  abline( h = seq(60,140,20), lty = 3, col = "gray45")
  if(mm == length(mods)) legend("bottom", ncol = 5, bty = "n", col = cols[1:5], legend = periods, lwd = lwd)
}
dev.off()

#### all models together
#cols <- c("#FF4500", "#FFD700", "#FF69B4", "#4682B4", "#2E8B57")
png(paste0(recon_dir,"/figures/mcomp_allmods_20yrperiods_diffs_",Sys.Date(),".png"), width = 1500, res = 300, height = 2000, pointsize = 7)
par(mfrow = c(length(mods),1), mar = c(0.5,2,0.5,1), oma= c(3,1,1,1))
for (mm in 1:length(mods)){
  plot(1763:2020,doymat[,1,1], type = "l", ylim = c(-10, 10), lwd = lwd, col = "white", xaxt = ifelse(mm == 7, "s","n"))
  for(i in 1:length(periods)) points(1763:2020,doymat[,i,mm] - doymat[,5,mm], type = "l", lwd = lwd, col = cols[i])
  mtext(side = 3, mods[mm], line = -1.2)
  abline( h = seq(-20,20,5), lty = 3, col = "gray45", lwd = 0.5)
  if(mm == length(mods)) legend("bottom", ncol = 5, bty = "n", col = cols[1:5], legend = periods, lwd = lwd)
}
dev.off()

#### Boxplot plots
doymat_long <- data.frame(matrix(doymat, ncol = prod(dim(doymat)[2:3])))
colnames(doymat_long) <- paste0(rep(mods, each = length(periods)), 1:length(periods))

doymat_long <- data.frame(model = rep(mods,each = 258 * length(periods)),
                          per = rep(rep(1:length(periods),each = 258),length(mods)),
                          year = rep(1763:2020,length(periods)*length(mods)),
                          doy = as.vector(doymat))

png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_20yrsperiods_18501950_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(doymat_long[doymat_long$year %in% 1850:1950,], aes(x=model, y=doy, fill = factor(per)) )+ # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = periods, name = "Calibration periods") +
  theme_bw() +
  theme(plot.title = element_text(size=9),
        legend.title = element_text(size=9)) +
  ggtitle("Distribution of model and period 1850 to 1950") +
  xlab("") + ylab("")
dev.off()


png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_20yrsperiods_17641950_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(doymat_long[doymat_long$year %in% 1764:1950,], aes(x=model, y=doy, fill = factor(per)) )+ # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = periods, name = "Calibration periods") +
  theme_bw() +
  theme(plot.title = element_text(size=9),
        legend.title = element_text(size=9)) +
  ggtitle("Distribution of model and period 1764 to 1950") +
  xlab("") + ylab("")
dev.off()

png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_20yrsperiods_18001900_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(doymat_long[doymat_long$year %in% 1800:1900,], aes(x=model, y=doy, fill = factor(per)) )+ # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = periods, name = "Calibration periods") +
  theme_bw() +
  theme(plot.title = element_text(size=9),
        legend.title = element_text(size=9)) +
  ggtitle("Distribution of model and period 1800 to 1900") +
  xlab("") + ylab("")
dev.off()


png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_20yrsperiods_19912020_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(doymat_long[doymat_long$year %in% 1991:2020,], aes(x=model, y=doy, fill = factor(per)) )+ # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = F,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = periods, name = "Calibration periods") +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    legend.title = element_text(size=9)
  ) +
  ggtitle("Distribution of model and period 1991 to 2020") +
  xlab("") + ylab("")
dev.off()

####### compare to observations
predicted_values <- lapply(mod_outs,function(x) {
  dat <- data.frame(sapply(1:length(x$modelled), function(mm) x$modelled[[mm]]$predicted_values))
  colnames(dat) <- mods
  return(dat)
})

obs_list <- lapply(stn_lists, function(x){data.frame(site = x$site, yr= x$year, doy =x$transition_dates)})
obs_values <- dplyr::bind_rows(obs_list, .id = "periods")
#measured_values <- unlist(lapply(mod_outs,function(x) {x$measured}))
pred_tab <- dplyr::bind_rows(predicted_values, .id = "periods")
pred_tab <- melt(cbind(pred_tab,obs_values$doy))
pred_tab$periods <- paste0("P",stringr::str_sub(pred_tab$periods,-8,-1))

png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_calibrationperiods_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(pred_tab, aes(x=periods , y=value, fill= factor(variable))) + # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = c(mods,"measured"), name = "Models") +
  theme_bw() +
  theme(
    plot.title = element_text(size=8),
    legend.title = element_text(size=8)
  ) +
  ggtitle("Model prediction in calibration period") +
  xlab("") + ylab("day-of-year")
dev.off()

### calculate some errors
test_out <- lapply(1:5,function(per){
  out <- sapply(1:ncol(predicted_values[[per]]), function(x){
  predicted_values[[per]][,x] - obs_list[[per]]$doy
    })
  rownames(out) <- paste0(obs_list[[per]]$site,"-",obs_list[[per]]$yr)
  return(out)
}
)


png(paste0(recon_dir,"/figures/mcomp_allmods_boxplots_calibrationperiods_",Sys.Date(),".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
ggplot(pred_tab, aes(x=periods , y=value, fill= factor(variable))) + # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T, labels = c(mods,"measured"), name = "Models") +
  theme_bw() +
  theme(
    plot.title = element_text(size=8),
    legend.title = element_text(size=8)
  ) +
  ggtitle("Model prediction in calibration period") +
  xlab("") + ylab("day-of-year")
dev.off()
