# Generate model drivers from climate reconstruction data

# load packages
library(phenor)

# select phenophases for calibration
classes <- "1234_longseries_filled"
pv <- "mpica13d"

# where to save the data
output_dir <- "data/pheno_net/"

# create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Generating complete driver data (all years, no subsets)")
# --- run calibration for alternating years ----
message(paste0("-- processing: ", pv))

# set filenames for driver data
meteo_file <- file.path("data-raw/pheno_prep/",
                        paste0("meteo_dat_",pv,"_class",classes,".rds")
)
phenology_file <- file.path("data-raw/pheno_prep/",
                            paste0("pheno_dat_",pv,"_class",classes,".rds")
)

# read pheno data and meteorological data
meteo_dat <- readRDS(meteo_file)
pheno_dat <- readRDS(phenology_file)

# remove nas
naind <- is.na(pheno_dat$doy)
pheno_dat <- pheno_dat[!naind,]
meteo_dat <- meteo_dat[!naind,]

# split out station names
stns_pv <- unique(pheno_dat$nat_abbr)

# create empty list to hold all data
stn_list <- list()

# convert data to phenor format
for(stn in 1:length(stns_pv)){

  # get id to process
  id <- stns_pv[stn]

  ind <- which(pheno_dat$nat_abbr == id & pheno_dat$reference_year < 2021)
  if(length(ind) > 1)ltm <- apply(meteo_dat[rownames(meteo_dat) == id,], 2, mean) else {ltm <- meteo_dat[rownames(meteo_dat) == id,]}

  # create meteo data subset
  meteo_subset <- as.matrix(meteo_dat[rownames(meteo_dat) == id,])
  if(dim(meteo_subset)[1] != 365) meteo_subset <- t(meteo_subset)
  # check if this is a matrix
  # if not subsetting the array
  # failed and the site does not exist
  # (and should be skipped)
  if(!length(meteo_subset) > 0){
    message(paste0("--- missing site: ", stn))
    stn_list[[stn]] <- NULL
  } else {

    locs <- as.numeric(unlist(pheno_dat[ind[1],c("obsLat","obsLon")]))

    stn_list[[stn]]  <- list(
      site = id,
      location = locs,
      doy = c(-102:-1,1:263),
      ltm = ltm, # long term mean temperature for a given location / ignored
      transition_dates = pheno_dat$doy[ind],
      year = pheno_dat$reference_year[ind],
      Ti = meteo_subset,
      Tmini = NULL, # tmin and max ignored / allowed to be empty
      Tmaxi = NULL, #
      Li = matrix(
        rep(daylength(c(264:365,1:263),locs[1]),
            length(ind)),
        nrow = 365,
        ncol = length(ind)
      )
    )
  }
}

# add names to the list elements
names(stn_list) <- stns_pv

# flatten the phenocam nested list format
stn_list <- pr_flatten(stn_list)

# construct filename
filename <- file.path(output_dir, paste0(pv, "_full_model_drivers_class",classes,".rds"))

# save data as a compressed RDS file
saveRDS(
  stn_list,
  filename,
  compress = "xz"
)

