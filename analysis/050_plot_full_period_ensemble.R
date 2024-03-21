# plot the recons based on model an periods

# packages
library(ggplot2)
library(dplyr)
library(viridis)
library(reshape2)

# set variables
classes <- "1234_longseries_filled"
burnin <- 4000
iter <- 36000
date <- "2024-03-11"
pv <- "mprua65d"

# directories
base_dir <- file.path("data/calibrations",paste0("bt_run_",date,"_class",classes,"_",burnin,"burnin_",iter,"iter"))
recon_dir <- file.path(base_dir,pv)

# read model data
files <- list.files(file.path(base_dir,pv), pattern = "model_comparison", full.names = T)
mod_nams <- paste0("P",stringr::str_sub(files,-12,-5))
mod_outs <- lapply(files, readRDS)
names(mod_outs) <- mod_nams

# read observation data by period
stn_files <- list.files(file.path(base_dir,pv), pattern = "stn_list_train", full.names = T)
stn_lists <- lapply(stn_files, readRDS)
names(stn_lists) <- mod_nams

# read historical phenological predictions
doymat <- readRDS(paste0(recon_dir,"/liestal_estimate_",pv,"_periods_model_ensemble_",classes,".rds"))
dd <- dim(doymat)

mods <-  names(mod_outs$P19611980$modelled)
periods <- paste0(stringr::str_sub(mod_nams,2,5),"-",stringr::str_sub(mod_nams,6,10))

cols <- c("#004177", "#005694", "#0075bc", "#008cda", "#00a5ff",  "#660086", "#9200ad", "#b600d0","#9932CC", "#e6e6ff")
lwd <- 1

# 1764-1960 time series of all models together
png(paste0(recon_dir,"/figures/mcomp_allmods_20yrperiods.png"), width = 1500, res = 300, height = 2000, pointsize = 7)
par(mfrow = c(length(mods),1), mar = c(0.5,2,0.5,1), oma= c(3,1,1,1))
for (mm in 1:length(mods)){
  plot(1763:1960,doymat[,1,1,1], type = "l", ylim = c(80, 135), lwd = lwd, col = "white", xaxt = ifelse(mm == 7, "s","n"))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
  for(i in 1:length(periods))
    sapply(1:10, function(e) {points(1763:1960,doymat[,e,i,mm], type = "l", lwd = lwd, col = cols[i])})
  mtext(side = 3, mods[mm], line = -1.3)
  abline( h = seq(60,140,20), lty = 3, col = "gray45")
  if(mm == length(mods)) legend("bottom", ncol = 5, bty = "n", col = cols[1:length(periods)], legend = periods, lwd = lwd + 0.1)
}
dev.off()

# 1764-1960 difference time series of all models together
png(paste0(recon_dir,"/figures/mcomp_allmods_20yrperiods_diffs.png"), width = 1500, res = 300, height = 2000, pointsize = 7)
par(mfrow = c(length(mods),1), mar = c(0.5,2,0.5,1), oma= c(3,1,1,1))
for (mm in 1:length(mods)){
  plot(1763:1960,doymat[,1,1,1], type = "l", ylim = c(-15, 15), lwd = lwd, col = "white", xaxt = ifelse(mm == 7, "s","n"))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
  for(i in 1:(length(periods)-1))
    sapply(1:10, function(e) {points(1763:1960,doymat[,e,i,mm] - doymat[,e,length(periods),mm], type = "l", lwd = lwd, col = cols[i])})
  mtext(side = 3, mods[mm], line = -1.3)
  abline( h = seq(-20,20,5), lty = 3, col = "gray45", lwd = 0.5)
  if(mm == length(mods)) legend("bottom", ncol = 5, bty = "n", col = cols[1:length(periods)], legend = periods, lwd = lwd + 0.1)
  if(mm == 1) mtext(paste0("Difference to period ",periods[length(periods)]),3, )
}
dev.off()

# Boxplots
doymat_long <- data.frame(model = rep(mods,each = prod(dd[c(1,2,3)])),
                          per = rep(rep(1:length(periods),each = dd[1] * dd[2]),length(mods)),
                          year = rep(1763:1960,length(periods)*length(mods)*dd[2]),
                          ens = rep(rep(1:10,each = dd[1]),prod(dd[3:4])),
                          doy = as.vector(doymat))

aggregate(doymat_long$doy, list(doymat_long$per,doymat_long$model), median, na.rm = T)
level_periods <- c(1,4,6,8,10,2,5,7,9,3)

# not the change of colors with respct to plots above
png(paste0(recon_dir,"/figures/mcomp_allmods_",pv,"_boxplots_allperiods_17641950.png"), width = 2500, res = 300, height = 1000, pointsize = 5)
ggplot(doymat_long[doymat_long$year %in% 1764:1950,], aes(x=model, y=doy, fill = factor(per, level = level_periods)) )+ # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.2) +
  scale_fill_manual(values = cols, labels = periods[level_periods], name = "Calibration periods") +
  theme_bw() +
  theme(plot.title = element_text(size=8),
        legend.title = element_text(size=8),
        axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +
  ggtitle("Distribution of model and period 1764 to 1950") +
  xlab("") + ylab("")
dev.off()

# compare to observations
predicted_values <- lapply(mod_outs,function(x) {
  dat <- data.frame(sapply(1:length(x$modelled), function(mm) x$modelled[[mm]]$predicted_values))
  colnames(dat) <- mods
  return(dat)
})

obs_list <- lapply(stn_lists, function(x){data.frame(site = x$site, yr= x$year, doy =x$transition_dates)})
obs_values <- dplyr::bind_rows(obs_list, .id = "periods")
pred_tab <- dplyr::bind_rows(predicted_values, .id = "periods")
pred_tab <- melt(cbind(pred_tab,obs_values$doy))
pred_tab$periods <- paste0("P",stringr::str_sub(pred_tab$periods,-8,-1))

png(paste0(recon_dir,"/figures/mcomp_allmods_",pv,"_boxplots_calibrationperiods.png"), width = 2400, res = 300, height = 1000, pointsize = 5)
ggplot(pred_tab, aes(x=periods , y=value, fill= factor(variable))) + # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.2) +
  scale_fill_viridis(discrete= T, labels = c(mods,"measured"), name = "Models") +
  theme_bw() +
  theme(
    plot.title = element_text(size=8),
    legend.title = element_text(size=8),
    axis.text=element_text(size=7),
    axis.title=element_text(size=8)
  ) +
  ggtitle("Model predictions in each calibration period") +
  xlab("") + ylab("day-of-year")
dev.off()

# calculate differences
diff_out <- lapply(1:length(periods),function(per){
  out <- data.frame(sapply(1:ncol(predicted_values[[per]]), function(x){
    predicted_values[[per]][,x] - obs_list[[per]]$doy
  }))
  colnames(out) <- mods
  return(out)
})

names(diff_out) <- names(predicted_values)

diff_tab <- dplyr::bind_rows(diff_out, .id = "periods")
diff_tab <- melt(diff_tab)
diff_tab$periods <- paste0("P",stringr::str_sub(diff_tab$periods,-8,-1))

png(paste0(recon_dir,"/figures/mcomp_allmods_",pv,"_boxplots_calibrationperiods_differences.png"), width = 2400, res = 300, height = 1000, pointsize = 5)
ggplot(diff_tab, aes(x=periods , y=value, fill= factor(variable))) + # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.2) +
  scale_fill_viridis(discrete= T, labels = c(mods,"measured"), name = "Models") +
  theme_bw() +
  theme(
    plot.title = element_text(size=8),
    legend.title = element_text(size=8),
    axis.text=element_text(size=7),
    axis.title=element_text(size=8)
  ) +
  ggtitle("Model prediction in calibration period") +
  xlab("") + ylab("day-of-year")
dev.off()
