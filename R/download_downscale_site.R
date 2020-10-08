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
                                    downscale,
                                    overwrite,
                                    noaa_var_names = NULL,
                                    model_name,
                                    model_name_ds,
                                    output_directory, 
                                    ensemble=10, 
                                    tend=40 #10 days 4 sim per day
                                    ) {
  
  model_dir <- file.path(output_directory, model_name)


  #Create dimensions of the noaa forecast
  lon.dom <- seq(0, 359, by = 0.5) #domain of longitudes in model (1 degree resolution)
  lat.dom <- seq(-90, 90, by = 0.5) #domain of latitudes in model (1 degree resolution)

  #Convert negetive longitudes to degrees east
  if(lon_list[site_index] < 0){
    lon_east <- 360 + lon_list[site_index]
  }else{
    lon_east <- lon_list[site_index]
  }

  urls.out <- tryCatch(rNOMADS::GetDODSDates(abbrev = "gens_bc"),
                       error = function(e){
                         warning(paste(e$message, "NOAA Server not responsive"),
                                 call. = FALSE)
                         return(NA)
                       },
                       finally = NULL)

  #urls.out <- list()

  urls.out$model <- "gefs"
  urls.out$date <- urls.out$date[which(lubridate::as_date(urls.out$date) >= lubridate::as_date("2020-09-24"))]
  urls.out$url <- paste0("https://nomads.ncep.noaa.gov:443/dods/gefs/gefs",urls.out$date)


  if(is.na(urls.out[1])) stop()

  urls.out$date_formated <- lubridate::as_date(urls.out$date)

  if(forecast_date == "latest"){
    url_index <- length(urls.out$url)
    previous_day_index <- url_index - 1
  }else if(!forecast_date == "all"){
    url_index <- which(urls.out$date_formated == lubridate::as_date(forecast_date))
  }else{
    url_index <- 1:length(urls.out$url)
  }

  for(i in url_index){

    model.url <- urls.out$url[i]
    start_date <- urls.out$date[i]

    model_list <- c("gefs_pgrb2ap5_all_00z", "gefs_pgrb2ap5_all_06z", "gefs_pgrb2ap5_all_12z", "gefs_pgrb2ap5_all_18z")
    model_hr <- c(0, 6, 12, 18)

    model.runs <- tryCatch(rNOMADS::GetDODSModelRuns(model.url),
                           error = function(e){
                             warning(paste(e$message, "skipping", model.url),
                                     call. = FALSE)
                             return(NA)
                           },
                           finally = NULL)
    if(is.na(model.runs)[1]) next

    avail_runs <- model.runs$model.run[which(model.runs$model.run %in% model_list)]
    avail_runs_index <- which(model_list %in% avail_runs)

    model_list <- model_list[avail_runs_index]
    model_hr <- model_hr[avail_runs_index]
    if(forecast_time != "all" & forecast_time != "latest"){
      hour_index <- which(model_hr %in% forecast_time)
      model_list <- model_list[hour_index]
    }else if(forecast_time == "latest"){
      hour_index <- which.max(model_hr)
      model_list <- model_list[hour_index]
    }


    for(m in 1:length(model_list)){

      run_hour <- stringr::str_sub(model_list[m], start = 19, end = 20)
      start_time <- lubridate::as_datetime(start_date) + lubridate::hours(as.numeric(run_hour))
      end_time <- start_time + lubridate::days(16)

      model_site_date_hour_dir <- file.path(model_dir, site_list[site_index], start_date,run_hour)
      if(!dir.exists(model_site_date_hour_dir)){
        dir.create(model_site_date_hour_dir, recursive=TRUE, showWarnings = FALSE)
      }

      identifier <- paste(model_name, site_list[site_index], format(start_time, "%Y-%m-%dT%H"),
                          format(end_time, "%Y-%m-%dT%H"), sep="_")

      #Check if already downloaded
      if(length(list.files(model_site_date_hour_dir)) != 31){

        print(paste("Downloading", site_list[site_index], format(start_time, "%Y-%m-%dT%H")))

        if(is.na(model.runs)[1]) next

        #check if available at NOAA
        if(model_list[m] %in% model.runs$model.run){

          model.run <- model.runs$model.run[which(model.runs$model.run == model_list[m])]

          cf_var_units1 <- c("K", "Pa", "1", "Wm-2", "Wm-2", "kgm-2s-1", "1", "1", "ms-1")  #Negative numbers indicate negative exponents

          noaa_data <- list()

          download_issues <- FALSE
  
          for(j in 1:length(noaa_var_names)){

            #For some reason rNOMADS::GetDODSDates doesn't return "gens" even
            #though it is there
            curr_model.url <- model.url

            lon <- which.min(abs(lon.dom - lon_east)) - 1 #NOMADS indexes start at 0
            lat <- which.min(abs(lat.dom - lat_list[site_index])) - 1 #NOMADS indexes start at 0

 
            noaa_data[[j]] <- tryCatch(rNOMADS::DODSGrab(model.url = curr_model.url,
                                                         model.run = model.run,
                                                         variables	= noaa_var_names[j],
                                                         time = c(0, tend),
                                                         lon = lon,
                                                         lat = lat,
                                                         ensembles=c(0, ensemble)),
                                       error = function(e){
                                         warning(paste(e$message, "skipping", curr_model.url, model.run, noaa_var_names[j]),
                                                 call. = FALSE)
                                         return(NA)
                                       },
                                       finally = NULL)

            if(is.na(noaa_data[[j]][1])){
              download_issues <- TRUE
            }else if(length(unique(noaa_data[[j]]$value)) == 1){
              #ONLY HAS 9.999e+20 which is the missing value
              download_issues <- TRUE
            }


            #For some reason it defaults to the computer's time zone, convert to UTC
            noaa_data[[j]]$forecast.date <- lubridate::with_tz(noaa_data[[j]]$forecast.date,
                                                               tzone = "UTC")
          }

          if(download_issues == TRUE){
            warning(paste("Error downloading one of the variables: ", curr_model.url, model.run))
            next
          }

          names(noaa_data) <- noaa_var_names

          
          noaa_df <- names(noaa_data) %>%
            map_dfc(~ data.frame(noaa_data[[.x]]$value) %>%
                      `colnames<-`(c(.x))
            ) %>%
            mutate(ensmbles=noaa_data[[1]]$ensembles,
                   time=noaa_data[[1]]$forecast.date)
          
          if(all(c("ugrd10m","vgrd10m") %in% names(noaa_df))){
            noaa_df <- noaa_df %>%
              mutate(wind_speed = sqrt(noaa_df$ugrd10m^2 + noaa_df$vgrd10m^2)) %>%
              dplyr::select(-ugrd10m, -vgrd10m)
          }
          
          if("dlwrfsfc" %in% names(noaa_df)){
            noaa_df <- noaa_df %>%
              mutate(dlwrfsfc = replace(dlwrfsfc, dlwrfsfc == 9.999e+20, NA))
          }
          
          if("dswrfsfc" %in% names(noaa_df)){
            noaa_df <- noaa_df %>%
              mutate(dswrfsfc = replace(dswrfsfc, dswrfsfc == 9.999e+20, NA))
          }
          if("apcpsfc" %in% names(noaa_df)){
            noaa_df <- noaa_df %>%
              mutate(apcpsfc = replace(apcpsfc, apcpsfc == 9.999e+20, NA))
          }
          
          


      return(noaa_df)
        }
      }else{
        print(paste("Existing", site_list[site_index], start_time))
      }
    }
  }
}
