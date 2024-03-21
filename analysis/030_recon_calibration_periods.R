# Reconstruct each calibration periods based on all other calibration periods

# packages
library(phenor)
library(reshape2)
library(ggplot2)
library(viridis)
library(dplyr)

# set variables
pv <- "mprua65d"
classes <- "1234_longseries_filled"
burnin <- 4000
iter <- 36000
date <- "2024-03-11"
selpers <- c("19611980|19711990|19812000|19912010|20012020|19612020")

# directories
base_dir <- file.path(
  "data/calibrations",
  paste0("bt_run_",date,"_class",classes,"_",burnin,"burnin_",iter,"iter")
)

recon_dir <- file.path(base_dir,pv)

# read model comparison files
files <- list.files(file.path(base_dir,pv), pattern = "model_comparison", full.names = T)
files <- files[grepl(pattern = selpers, files)]
cal_nams <- paste0("P",stringr::str_sub(files,-12,-5))
mod_outs <- lapply(files, readRDS)
names(mod_outs) <- cal_nams

# read observation files
stn_files <- list.files(file.path(base_dir,pv), pattern = "stn_list_train", full.names = T)
stn_files <- stn_files[grepl(pattern = selpers,stn_files)]
stn_lists <- lapply(stn_files, readRDS)
names(stn_lists) <- cal_nams
per_nams <- paste0(stringr::str_sub(cal_nams,2,5),"-",stringr::str_sub(cal_nams,6,10))

per_stn_list <- lapply(1:length(mod_outs), function(per){

  per_help <- mod_outs[[per]]
  mods <- names(per_help$modelled)

  stn_lists_out <- lapply(stn_lists, function(ss){

    doyout <- sapply(1:length(per_help$modelled), function(mm){
      # get model and parameters
      opt_params <- per_help$modelled[[mm]]$parameters
      func <- get(mods[mm])
      doyvec <- func(par = opt_params, data = ss)
    })

    measured_values <- ss$transition_dates
    doyout <- data.frame(cbind(doyout, measured_values))
    colnames(doyout) <- c(mods,"obs")
    return(doyout)
  })

  stn_lists_tab <- bind_rows(stn_lists_out, .id = "periods")
  stn_lists_tab <- melt(stn_lists_tab)
  stn_lists_tab[stn_lists_tab$value > 365,"value"] <- NA
  stn_lists_tab$periods[stn_lists_tab$periods==cal_nams[per]] <- paste0("calp_",cal_nams[per])

  return(stn_lists_tab)

})

names(per_stn_list) <- per_nams
mods <- names(mod_outs$P19611980$modelled)

# plot for each period model predictions as boxplots
for(per in 1:length(per_stn_list)){
  png(paste0(recon_dir,"/figures/mcomp_",pv,"_boxplots_calibrated_in_period_",per_nams[per],".png"), width = 2000, res = 300, height = 1000, pointsize = 5)
  p <- ggplot(per_stn_list[[per]], aes(x=periods , y=value, fill= factor(variable))) + # fill=name allow to automatically dedicate a color for each group
    stat_boxplot(geom ='errorbar',linewidth = 0.25) +
    geom_boxplot(notch = T,linewidth = 0.25,outlier.size=0.5) +
    scale_fill_viridis(discrete= T, labels = c(mods,"measured"), name = "Models") +
    theme_bw() +
    theme(
      plot.title = element_text(size=8),
      legend.title = element_text(size=8),
      axis.text=element_text(size=7),
      axis.title=element_text(size=8)
    ) +
    ggtitle(paste0("Model prediction of all calibration periods based models calibrated in the period ", per_nams[per])) +
    xlab("") + ylab("day-of-year")
  print(p)
  dev.off()
}

