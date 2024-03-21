# Calibration routine using different periods
# over the whole dataset

try(detach("package:phenor", unload = TRUE))
# load packages
library(phenor)
library(dplyr)
library(BayesianTools)

# select phenophases for calibration
pv <- "mfags13d"
classes <- "1234_longseries_filled"
burnin <- 4000
iter <- 36000

# create directory for new calibration with simulated annealing
base_dir <- file.path(
  "data/calibrations/",
  paste0("bt_run_",Sys.Date(),"_class",classes,"_",burnin,"burnin_",iter,"iter")
)

dir.create(
  base_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# select seeds, model, years
mods <-  c("LIN","TT", "TTs", "PTT", "PTTs","M1","M1s","AT")

# list driver files
driver_files <- list.files("data/pheno_net/","*.rds",full.names = TRUE)

message("Calibrating models:")

# read in the full driver files
drivers <- readRDS(driver_files[grepl(paste0(pv,".*",classes,".rds"), driver_files)])

# set periods
periods <- list(c(1961,1980),c(1971,1990),c(1981,2000),c(1991,2010),c(2001,2020),c(1961,1990),c(1971,2000),c(1981,2010),c(1991,2020),c(1961,2020))

# --- run calibration for alternating years ----
lapply(periods, function(per) {

  message(paste0("-- calibrating: ", pv, per[1], "-", per[2]))

  # create selection criteria
  # for training
  selection <- data.frame(
    site = drivers$site,
    year = drivers$year
  ) |>
    mutate(
      train = (year %in% per[1]:per[2])
    )

  # select training years
  selection <- selection |>
    dplyr::group_by(site) |>
    mutate(
      train = ifelse(
        length(which(train)) > 1 & train == TRUE,
        TRUE,
        FALSE
      )
    )

  # subset training and testing datasets
  stn_list_train <- pr_fm_subset(drivers, selection$train)

  print("fit models")

  mod_out <- pr_fit_comparison(
    random_seeds = 1,
    models = mods,
    data = stn_list_train,
    method = "bayesiantools",
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = burnin,
        iterations = iter
      )
    )
  )

  geldiag <- lapply(mod_out$modelled, function(x){gelmanDiagnostics(x$opt_out, plot = F)})


  ### set 9999 to NA for evaluation
  for(ii in 1:length(mod_out$modelled)){
    mod_out$modelled[[ii]]$predicted_values[mod_out$modelled[[ii]]$predicted_values > 1000] <- NA
  }

  # create output directory
  dir.create(
    file.path(base_dir, pv),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # save data as compressed RDS files
  saveRDS(
    mod_out,
    file = file.path(base_dir,pv,paste0("model_comparison_",pv,"_",per[1],per[2],".rds")),
    compress = "xz"
  )


  saveRDS(
    geldiag,
    file = file.path(base_dir,pv,paste0("gelman_diagnostics_",pv,"_",per[1],per[2],".rds")),
    compress = "xz"
  )

  saveRDS(
    stn_list_train,
    file = file.path(base_dir,pv,paste0("stn_list_train_",pv,"_",per[1],per[2],".rds")),
    compress = "xz"
  )

  # create directory for plot
  dir.create(
    file.path(base_dir, pv, "figures"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # plot the model convergence for every model and period
  for(mm in 1:length(mods)){
    pdf(file.path(base_dir, pv, "figures",paste0("mod_covergence_",mods[mm],"_",per[1],per[2],".pdf")))
    plot(mod_out$modelled[[mm]]$opt_out)
    dev.off()
  }

})
