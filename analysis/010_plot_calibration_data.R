# plot climatologies of calibration data based on different periods
rm(list=ls())

library(dplyr)
library(phenor)

# library(RColorBrewer)
library(ggplot2)
library(viridis)

classes <- "1234_longseries_filled"
pv <- "mlard13d"
periods <- list(c(1961,1980),c(1971,1990),c(1981,2000),c(1991,2010),c(2001,2020))

# create directory for plots by phenophase
plot_dir <- file.path("analysis/figures/pheno_clim",pv)
dir.create(plot_dir,showWarnings = FALSE,recursive = TRUE)

# read driver files
driver_files <- list.files("data/pheno_net/","*.rds", full.names = TRUE)
drivers <- readRDS(driver_files[grepl(paste0(pv,".*",classes,".rds"), driver_files)])

stn_lists <- lapply(periods, function(per){

  selection <- data.frame(
    site = drivers$site,
    year = drivers$year
  ) |>
    mutate(
      train = (year %in% per[1]:per[2])
    )

  drivers_per <- pr_fm_subset(drivers, selection$train)
  df <- data.frame(site = drivers_per$site,
                   date= drivers_per$transition_dates,
                   yr = drivers_per$year,
                   lat = drivers_per$location[1,],
                   lon = drivers_per$location[2,])
  return(df)
})

pernam <- sapply(periods,function(x)paste0("P",x[1],x[2]))
names(stn_lists) <- pernam

stn_lists_tab <-  bind_rows(stn_lists, .id = "periods")

png(paste0(plot_dir,"/",pv,"_phenophase_distribution_per_period_boxplot_",classes,"_",Sys.Date(),".png"),
    width = 1100, height = 600, res = 300, pointsize=3)
ggplot(stn_lists_tab, aes(x=periods , y=date, fill = periods)) + # fill=name allow to automatically dedicate a color for each group
  stat_boxplot(geom ='errorbar',linewidth = 0.25) +
  geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
  scale_fill_viridis(discrete= T,guide="none") +
  theme_bw() +
  theme(
    plot.title = element_text(size=8),
    legend.title = element_text(size=8),
    axis.title.y=element_text(size=8)
  ) +
  ggtitle("Phenophase distribution per period") +
  xlab("") + ylab("day of year")
dev.off()

png(paste0(plot_dir,"/",pv,"_timeseries_per_period_",classes,"_",Sys.Date(),".png"), width = 1100, height = 1700, res = 300, pointsize=7)
par(mfrow= c(length(periods),1), mar= c(2,2,1,1))
lapply(stn_lists, function(per){
  print(sites <- unique(per$site))
  cols <- viridis(n = length(sites))
  plot(per$yr, per$date, col = "white", ylim = c(75,155), xlab = "")
  abline(h = seq(10,150,10), lty = 3, col = "gray45")
  for(ss in 1:length(sites)){
    points(per$yr[per$site==sites[ss]], per$date[per$site==sites[ss]], type = "p", col = cols[ss], pch =19)
  }
  climyr <- aggregate(per$date,list(per$yr),mean)
  points(climyr$Group.1,climyr$x, type = "l")
  points(climyr$Group.1,climyr$x, type = "p",19)
})
dev.off()


brks <- seq(75,140,5)
cols <- rev(viridis(n = length(brks)-1))

png(paste0(plot_dir,"/",pv,"_mean_phenophase_per_period_map_",classes,"_",Sys.Date(),".png"), width = 1500, height = 600, res = 300, pointsize=5)
layout(matrix(1:6,nrow=2, byrow = T))
lons <- seq(5.7,10.5,0.1)
lats <- seq(45.5,48,0.1)
par(mar= c(1,1,1,1))
for(per in 1:length(periods)){
  image(lons,lats,matrix(NA,ncol = length(lats),nrow= length(lons)), asp = 1, bty = "n", xaxt = "n", yaxt = "n")
  maps::map("world",add= T)
  yrlen <- table(stn_lists[[per]]$site)
  clim <- aggregate(stn_lists[[per]][,2:5], list(stn_lists[[per]]$site), mean)
  print(range(clim$date))
  colp <- cols[cut(clim$date, brks)]
  points(clim$lon,clim$lat, pch = 19, col = colp, cex= yrlen/10 + 0.1)
  mtext(side=3, line = -1.4, pernam[per])
  if(per == length(periods)) {
    par(xpd = T)
    plot.new()
    legend("left",  ncol = 1,legend = levels(cut(clim$date, brks)), col = cols, pch = 19, bty = "n")
  }
}
dev.off()

# look at trends in these series
sites <- unique(drivers$site)

for(ss in 1:length(sites)){

  sind <- which(drivers$site == sites[ss])
  if(length(sind) > 8){
    site_dat <- data.frame(yr = drivers$year[sind], date = drivers$transition_dates[sind])
    site_dat<- site_dat[order(site_dat$yr),]

    png(paste0(plot_dir,"/trends_",pv,"_", sites[ss],".png"), width = 800, res = 200, pointsize = 6)
    plot(site_dat$yr, site_dat$date, type = "b", ylab = "doy", xlab = "year")
    lmmod <- lm(date~yr, data= site_dat)
    title(sites[ss])
    points(site_dat$yr,lmmod$fitted.values, type = "l")
    dev.off()
  }

}

# look at mean trend
drivers_test <- drivers
selind <- drivers_test$site == "AUR" & drivers_test$year < 1970
drivers_test$transition_dates[selind] <- NA
mean_site <- aggregate(drivers_test$transition_dates, list(drivers_test$year), mean, na.rm = T)
plot(mean_site$Group.1, mean_site$x, type = "l")
lmmod <- lm(x~Group.1, mean_site)
points(mean_site$Group.1, lmmod$fitted.values)
