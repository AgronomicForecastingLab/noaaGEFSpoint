##' @title Download and Downscale NOAA GEFS for a single site
##' @return None
##'
##' @param site_index, index of site_list, lat_list, lon_list to be downloaded
##' @param lat_list, vector of latitudes that correspond to site codes
##' @param lon_list, vector of longitudes that correspond to site codes
##' @param site_list, vector of site codes, used in directory and file name generation
##' @param downscale, logical specifying whether to downscale from 6-hr to 1-hr
##' @param overwrite, logical stating to overwrite any existing output_file
##' @param model_name, directory name for the 6-hr forecast, this will be used in directory and file name generation
##' @param model_name_ds, directory name for the 1-hr forecast, this will be used in directory and file name generation
##' @param output_directory, directory where the model output will be save
##' @export
##'
##' @author Quinn Thomas
##'
##'


download_downscale_site <- function(site_index,
                                    lat_list,
                                    lon_list,
                                    site_list,
                                    forecast_time = "all",
                                    forecast_date = "all",
                                    latest = FALSE,
                                    downscale,
                                    overwrite,
                                    model_name,
                                    model_name_ds,
                                    output_directory){

  model_dir <- file.path(output_directory, model_name)

  #Create dimensions of the noaa forecast
  lon.dom <- seq(0, 359, by = 1) #domain of longitudes in model (1 degree resolution)
  lat.dom <- seq(-90, 90, by = 1) #domain of latitudes in model (1 degree resolution)

  lon.dom.bc <- seq(0, 359, by = 0.5) #domain of longitudes in model (1 degree resolution)
  lat.dom.bc <- seq(-90, 90, by = 0.5) #domain of latitudes in model (1 degree resolution)

  #Convert negetive longitudes to degrees east
  if(lon_list[site_index] < 0){
    lon_east <- 360 + lon_list[site_index]
  }else{
    lon_east <- lon_list[site_index]
  }

  #For some reason rNOMADS::GetDODSDates doesn't return "gens" even
  #though it is there, use gens_bc to get the url but replace gens_bc with gens
  #below
  urls.out <- rNOMADS::GetDODSDates(abbrev = "gens_bc")

  if(latest){
    url_index <- length(length(urls.out$url))
  }else if(!forecast_date == "all"){
    url_index <- which(urls.out$date %in% forecast_date)
  }else{
  url_index <- 1:length(length(urls.out$url))
  }

  for(i in url_index){

    if(forecast_date == "all" | forecast_date %in% urls.out$date){

      model.url <- urls.out$url[i]
      start_date <- urls.out$date[i]

      model_list <- c("gep_all_00z", "gep_all_06z", "gep_all_12z", "gep_all_18z")
      model_hr <- c(0, 6, 12, 18)
      if(latest){
        model.runs <- GetDODSModelRuns(model.url)
        model_list <- model.runs$model.run[max(which(model.runs$model.run %in% model_list))]
      }else if(forecast_time != "all"){
        if(!forecast_time %in% c(0,6,12,18) & !latest){
          stop("forecast time not in avialable list c(0,6,12,18) in UTC")
        }
        model_list <- model_list[which(model_hr %in% forecast_time)]
      }

      model_site_date_dir <- file.path(model_dir, site_list[site_index], start_date)

      for(m in 1:length(model_list)){

        run_hour <- stringr::str_sub(model_list[m], start = 9, end = 10)
        start_time <- lubridate::as_datetime(start_date) + lubridate::hours(as.numeric(run_hour))
        end_time <- start_time + lubridate::days(16)

        model_site_date_hour_dir <- file.path(model_dir, site_list[site_index], start_date,run_hour)
        if(!dir.exists(model_site_date_hour_dir)){
          dir.create(model_site_date_hour_dir, recursive=TRUE, showWarnings = FALSE)
        }

        identifier <- paste(model_name, site_list[site_index], format(start_time, "%Y-%m-%dT%H"),
                            format(end_time, "%Y-%m-%dT%H"), sep="_")

        #Check if already downloaded
        if(length(list.files(model_site_date_hour_dir)) != 21){

          print(paste("Downloading", site_list[site_index], format(start_time, "%Y-%m-%dT%H")))

          model.runs <- rNOMADS::GetDODSModelRuns(model.url)

          #check if available at NOAA
          if(model_list[m] %in% model.runs$model.run){

            model.run <- model.runs$model.run[which(model.runs$model.run == model_list[m])]

            noaa_var_names <- c("tmp2m", "pressfc", "rh2m", "dlwrfsfc",
                                "dswrfsfc", "pratesfc",
                                "ugrd10m", "vgrd10m", "tcdcclm")

            noaa_var_names_model <- c("gens", "gens", "gens", "gens",
                                      "gens", "gens",
                                      "gens", "gens", "gens")

            #These are the cf standard names
            cf_var_names <- c("air_temperature", "air_pressure", "relative_humidity", "surface_downwelling_longwave_flux_in_air",
                              "surface_downwelling_shortwave_flux_in_air", "precipitation_flux", "eastward_wind", "northward_wind","cloud_area_fraction")

            #Replace "eastward_wind" and "northward_wind" with "wind_speed"
            cf_var_names1 <- c("air_temperature", "air_pressure", "relative_humidity", "surface_downwelling_longwave_flux_in_air",
                               "surface_downwelling_shortwave_flux_in_air", "precipitation_flux", "wind_speed","specific_humidity", "cloud_area_fraction")

            cf_var_units1 <- c("K", "Pa", "1", "Wm-2", "Wm-2", "kgm-2s-1", "ms-1", "1", "1")  #Negative numbers indicate negative exponents

            noaa_data <- list()

            for(j in 1:length(noaa_var_names)){

              #For some reason rNOMADS::GetDODSDates doesn't return "gens" even
              #though it is there
              curr_model.url <- stringr::str_replace(model.url, "gens_bc", noaa_var_names_model[j])

              if(noaa_var_names_model[j] == "gens_bc"){
                lon <- which.min(abs(lon.dom.bc - lon_east)) - 1 #NOMADS indexes start at 0
                lat <- which.min(abs(lat.dom.bc - lat_list[site_index])) - 1 #NOMADS indexes start at 0
              }else{
                lon <- which.min(abs(lon.dom - lon_east)) - 1 #NOMADS indexes start at 0
                lat <- which.min(abs(lat.dom - lat_list[site_index])) - 1 #NOMADS indexes start at 0
              }

              noaa_data[[j]] <- rNOMADS::DODSGrab(model.url = curr_model.url,
                                                  model.run = model.run,
                                                  variables	= noaa_var_names[j],
                                                  time = c(0, 64),
                                                  lon = lon,
                                                  lat = lat,
                                                  ensembles=c(0, 20))

              #For some reason it defaults to the computer's time zone, convert to UTC
              noaa_data[[j]]$forecast.date <- lubridate::with_tz(noaa_data[[j]]$forecast.date,
                                                                 tzone = "UTC")
            }

            names(noaa_data) <- cf_var_names

            forecast_noaa <- tibble::tibble(time = noaa_data$air_temperature$forecast.date,
                                            NOAA.member = noaa_data$air_temperature$ensembles,
                                            air_temperature = noaa_data$air_temperature$value,
                                            air_pressure= noaa_data$air_pressure$value,
                                            relative_humidity = noaa_data$relative_humidity$value,
                                            surface_downwelling_longwave_flux_in_air = noaa_data$surface_downwelling_longwave_flux_in_air$value,
                                            surface_downwelling_shortwave_flux_in_air = noaa_data$surface_downwelling_shortwave_flux_in_air$value,
                                            precipitation_flux = noaa_data$precipitation_flux$value,
                                            eastward_wind = noaa_data$eastward_wind$value,
                                            northward_wind = noaa_data$northward_wind$value,
                                            cloud_area_fraction = noaa_data$cloud_area_fraction$value)

            #9.999e+20 is the missing value so convert to NA
            forecast_noaa$surface_downwelling_longwave_flux_in_air[forecast_noaa$surface_downwelling_longwave_flux_in_air == 9.999e+20] <- NA
            forecast_noaa$surface_downwelling_shortwave_flux_in_air[forecast_noaa$surface_downwelling_shortwave_flux_in_air == 9.999e+20] <- NA
            forecast_noaa$precipitation_flux[forecast_noaa$precipitation_flux == 9.999e+20] <- NA
            forecast_noaa$cloud_area_fraction[forecast_noaa$cloud_area_fraction == 9.999e+20] <- NA

            forecast_noaa$specific_humidity <- rh2qair(rh = forecast_noaa$relative_humidity,
                                                       T = forecast_noaa$air_temperature,
                                                       press = forecast_noaa$air_pressure)

            #Calculate wind speed from east and north components
            forecast_noaa$wind_speed <- sqrt(forecast_noaa$eastward_wind^2 + forecast_noaa$northward_wind^2)

            #remove the east and north components
            forecast_noaa <- forecast_noaa %>% dplyr::select(-c("eastward_wind","northward_wind"))

            # Convert NOAA's total precipitation (kg m-2) to precipitation flux (kg m-2 s-1)
            #NOAA precipitation data is an accumulation over 6 hours.
            forecast_noaa$precipitation_flux <- udunits2::ud.convert(forecast_noaa$precipitation_flux,
                                                                     "kg m-2 hr-1",
                                                                     "kg m-2 6 s-1")  #There are 21600 seconds in 6 hours

            for (ens in 1:21) { # i is the ensemble number

              #Turn the ensemble number into a string
              if(ens< 10){
                ens_name <- paste0("0",ens)
              }else{
                ens_name <- ens
              }

              fname <- paste0(identifier,"_ens",ens_name,".nc")
              output_file <- file.path(model_site_date_hour_dir,fname)

              #Write netCDF
              write_noaa_gefs_netcdf(df = forecast_noaa,ens, lat = lat_list[site_index], lon_list[site_index], cf_units = cf_var_units1, output_file = output_file, overwrite = overwrite)

              if(downscale){
                #Downscale the forecast from 6hr to 1hr
                modelds_site_date_hour_dir <- file.path(output_directory,model_name_ds,site_list[site_index],start_date,run_hour)

                if(!dir.exists(modelds_site_date_hour_dir)){
                  dir.create(modelds_site_date_hour_dir, recursive=TRUE, showWarnings = FALSE)
                }

                identifier_ds <- paste(model_name_ds, site_list[site_index], format(start_time, "%Y-%m-%dT%H"),
                                       format(end_time, "%Y-%m-%dT%H"), sep="_")
                fname_ds <- file.path(modelds_site_date_hour_dir, paste0(identifier_ds,"_ens",ens_name,".nc"))

                #Run downscaling
                temporal_downscale(input_file = output_file, output_file = fname_ds, overwrite = overwrite)
              }
            }
          }
        }else{
          print(paste("Existing", site_list[site_index], start_time))
        }
      }
    }
  }
}