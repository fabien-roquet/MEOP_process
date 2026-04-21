import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import gsw
import cartopy.crs as ccrs


class CHLA_profile:
    """
    Functions for selecting, plotting and correcting from Quenching daytime/nighttime CHLA and light profiles.
    """
    def create_var(ds,varname,nprof):
        """
        Create a variable to plot profiles.
        :param ds: Dataset
        :param varname: variable of interest
        :param nprof: Number of profile 
        :return var: New variable in the dataset of nprof length 
        """
        var = xr.DataArray(ds[varname].isel(N_PROF=nprof), coords=[ds['PRES'].isel(N_PROF=nprof)], dims=['pres'])
        return var
    
    def select_200m(da):
        """
        Select the first 200m data of dataxarray
        :param da: Dataxarray
        :return: First 200m of the original dataxarray 
        """
        return da.where(ds.PRES<200,drop=True)

    def select_day(ds,day=True):
        """
        Select day profile (day=True) or night profiles (day=False)
        :param ds: Dataset
        :return: Dataset selection of day profile (day=True) or night profiles (day=False)
        """
        nprof = ds.is_day.where(ds.is_day==int(day),drop=True).N_PROF
        return ds.sel(N_PROF=nprof)

    def Chla_nozero(ds):
        """
        Select only non null CHLA profiles.
        :param ds: Dataset
        :return ds: Dataset with non null CHLA datarray values 
        """
        ds['CHLA_ADJUSTED'] = ds['CHLA_ADJUSTED'].where(ds['CHLA_ADJUSTED']!=0,np.nan)
        return ds

    def resample_day(ds):
        """
        Resample the dataset on a 24h mean.
        :param ds: Dataset
        :return: New resampled dataset 
        """
        ds2 = ds.rename({'N_PROF': 'DATE'})
        ds2['DATE'] = ds2.JULD
        return ds2.resample(DATE="24H").mean()

    def sunset_sunrise(time, lat, lon):
        """
        Calculates the local sunrise/sunset of the glider location.

        The function uses the Astral package to calculate the sunrise and sunset
        times using the date, latitude and longitude. The times are returned
        rather than day or night indices, as it is more flexible for the quenching
        correction.


        Parameters
        ----------
        time: numpy.ndarray or pandas.Series
            The date & time array in a numpy.datetime64 format.
        lat: numpy.ndarray or pandas.Series
            The latitude of the glider position.
        lon: numpy.ndarray or pandas.Series
            The longitude of the glider position.

        Returns
        -------
        sunrise: numpy.ndarray
            An array of the sunrise times.
        sunset: numpy.ndarray
            An array of the sunset times.

        """
        import astral as ast
        from astral.sun import sun
        from pandas import DataFrame

        df = DataFrame.from_dict(dict([("time", time), ("lat", lat), ("lon", lon)]))

        # set days as index
        df = df.set_index(df.time.values.astype("datetime64[D]"))

        # groupby days and find sunrise for unique days
        grp_avg = df.groupby(df.index).mean()
        date = grp_avg.index.to_pydatetime()

        sunrise_observer = []
        for i in range(len(grp_avg.lat)):
            sunrise_observer.append(
                ast.Observer(latitude=grp_avg.lat[i], longitude=grp_avg.lon[i])
            )

        sunrise, sunset = [], []
        for i in range(len(sunrise_observer)):
            sun_info = sun(sunrise_observer[i], date[i])
            sunrise.append(sun_info["sunrise"])
            sunset.append(sun_info["sunset"])

        grp_avg["sunrise"] = sunrise
        grp_avg["sunset"] = sunset
        # reindex days to original dataframe as night
        df_reidx = grp_avg.reindex(df.index).astype("datetime64[ns]")
        sunrise, sunset = df_reidx[["sunrise", "sunset"]].values.T

        return sunrise, sunset
    
    def is_day(time_of_day,time_of_sunrise,time_of_sunset):
        """
        Calculate a boolean variable of daytime for a given time of profile.
        :param time_of_day: Time of the profile (UTC datetime64 ?)
        :param time_of_sunrise: Time of sunrise of the day (UTC sunrise datetime64)
        :param time_of_sunset: Time of sunset of the day (UTC sunset datetime64)
        :return is_day: boolean variable daytime (is_day=True) or nighttime (is_day=False)
        """
        if time_of_sunrise<time_of_sunset:
            is_day = (time_of_day>time_of_sunrise) & (time_of_day<time_of_sunset)
        else:
            is_night = (time_of_day<time_of_sunrise) & (time_of_day>time_of_sunset)
            is_day = ~is_night
        return is_day

    def plot_ds(ds):
        """
        Plot CHLA_ADJUSTED, LIGHT_ADJUSTED, TEMP_ADJUSTED profiles of the dataset.
        :param ds: Dataset
        :return: Plot
        :raises: SunTimeException when there is no sunrise and sunset on given location and date
        """
        fig, axes = plt.subplots(ncols=1, nrows=3, sharey = True, figsize=(15,10))
        ds['CHLA_ADJUSTED'].dropna('N_PROF',how='all').T \
            .plot(yincrease=False,ax=axes[0],vmin=0,vmax=4)
        light_adjusted = (ds['LIGHT_ADJUSTED'].where(ds['LIGHT_ADJUSTED']!=0,np.nan)-ds['LIGHT_ADJUSTED'][:,18])
        light_adjusted.T.plot(yincrease=False,ax=axes[1],cmap=plt.cm.viridis,vmin=-4,vmax=0)
        ds['TEMP_ADJUSTED'].T.plot(yincrease=False,ax=axes[2],cmap=plt.cm.viridis,vmin=-2,vmax=10)
        plt.title(ds.smru_platform_code)

        for ax in axes:
            ax.set_ylim([0,400])
            ax.invert_yaxis()

    def plot_map(ds):
        """
        Plot a map of the tag trajectory around the Kerguelen plateau.
        :param ds: Dataset
        :return: Plot
        """
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.plot(ds.LONGITUDE,ds.LATITUDE,'.',transform=ccrs.PlateCarree(),label=fname)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
        ax.set_extent([30, 100, -60, -40], ccrs.PlateCarree())
        ax.stock_img()
        ax.coastlines()
        ax.gridlines(xlocs=ax.get_xticks(),ylocs=ax.get_yticks())

    def plot_prof_data(ds,nprof=200):
        """
        Plot profiles of temperature, salinity, chla and light.
        :param ds: Dataset
        :param nprof: Number of profiles to plot
        :return var: Plot
        """
        fig, axes = plt.subplots(ncols=4, nrows=1, sharey = True, figsize=(15,5))
        for n in range(nprof):
            create_var(ds,'TEMP_ADJUSTED',n).plot.line(y='pres', yincrease=False, ax=axes[0])
            create_var(ds,'PSAL_ADJUSTED',n).plot.line(y='pres', yincrease=False, ax=axes[1])
            create_var(ds,'CHLA_ADJUSTED',n).plot.line(y='pres', yincrease=False, ax=axes[2])
            create_var(ds,'LIGHT_ADJUSTED',n).plot.line(y='pres', yincrease=False, ax=axes[3])

    def to_day(ds,sunrise,sunset):
        """
        Convert dataset JULD data into datetime64 to create and add a is_day variable to the dataset.
        :param ds: Dataset
        :param sunrise: Time of sunrise of the day (UTC sunrise datetime64)
        :param sunset: Time of sunset of the day (UTC sunset datetime64)
        :return: Dataset with a is_day temporary variable
        """
        year = ds.JULD.dt.year.data
        month = ds.JULD.dt.month.data
        day = ds.JULD.dt.day.data
        hour = ds.JULD.dt.hour.data
        minute = ds.JULD.dt.minute.data
        second = ds.JULD.dt.second.data
        dates = [np.datetime64(datetime.datetime(year[i],month[i],day[i],hour[i],minute[i],second[i])) for i in range(ds.dims['N_PROF'])]

        for iprof in range(ds.dims['N_PROF']):
            if is_day(dates[iprof],sunrise[iprof],sunset[iprof]):
                ds['is_day'][iprof] = 1

    def mixed_layer_depth(ds, file_name):
        """
        Conpute the Mixed Layer Depth (MLD) as the depth where the density is higher than its value at 10 m by 0.02 kg m-3.
        :param ds: Dataset
        :param file_name: File name
        :return: Dataset with a MLD variable
        """
        N_PROF = ds.dims['N_PROF']
        mld = np.zeros(N_PROF)
        for i in range(N_PROF):
            pressure = ds.PRES[i,:].values
            temperature = ds.TEMP_ADJUSTED[i,:].values
            salinity = ds.PSAL_ADJUSTED[i,:].values
            index_t = np.isnan(temperature).argmin(axis=0)
            index_s = np.isnan(salinity).argmin(axis=0)
            temperature[:index_t]=temperature[index_t+1]
            salinity[:index_s]=salinity[index_s+1]
            if ~np.isnan(np.nansum(salinity)):
                density = gsw.sigma0(salinity,temperature)
                dens10 = density[9]
                ii = np.argmax(density-dens10>0.02)
                mld[i] = pressure[ii]
        mld[mld<5] = np.nan

        # save mld in netcdf file
        import netCDF4
        ncfile = netCDF4.Dataset(file_name,mode='a')
        if not('MLD' in ncfile.variables ):
            ncmld = ncfile.createVariable('MLD', np.float32, ('N_PROF',))
            ncmld.units = 'm'
            ncmld.long_name = 'mixed layer depth 0.02kg/m3'
            ncmld[:] = mld
        else:
            ncfile['MLD'][:] = mld
        ncfile.close()

    def photic_depth(par, nprof, depth, return_mask=False, ref_percentage=1):
        """
        Compute the euphotic depth according to the Glidertools function. 
        :param par: Light variable
        :param nprof: Number of profiles
        :param depth: Level variable
        :param ref_percentage: percentage of surface light
        :return: euphotic depth and attenuation coefficient
        """
        import numpy as np
        import pandas as pd
        from scipy.stats import linregress

        # Kd attenuation coefficient
        def dive_slope(par, depth):
            mask = ~(np.isnan(par) | np.isnan(depth))
            x, y = depth[mask], par[mask]
            slope = linregress(x, y).slope
            return slope

        # Percentage light depth
        def dive_light_depth(depth, slope):
            light_depth = np.exp((depth * -1) / (-1 / slope)) * 100.0
            ind = abs(light_depth - ref_percentage).argmin()
            euph_depth = depth[ind]
            return [euph_depth]

        slopes = []
        light_depths = []
        for d in range(len(nprof)):
            i = d
            zj = par[i,:]
            yj = depth[i,:]

            if (zj.all()==np.nan):
                slope = np.nan
            else:
                slope = dive_slope(zj, yj)
            light_depth = dive_light_depth(yj, slope)

            slopes += (slope,)
            light_depths += (light_depth,)

        slopes = pd.Series(slopes, index=nprof)
        light_depths = np.concatenate(light_depths)

        if not return_mask:
            light_depths = pd.Series(light_depths, index=nprof)

        return light_depths, slopes

    def quenching(ds):
        """
        Quenching correction according to Xing et al. 2018.
        :param ds: Dataset
        :return: dataset chlorophyll quenching corrected variable and create a temporay depth of light threshold variable
        """
        # Quenching according Xing et al 2018:
        ds['CHLA_CORRECTED'] = ds.CHLA_ADJUSTED.copy()

        iPAR = 1.17  # PAR Threshold (μmol quanta m-2 s-1)
        N_PROF = ds.dims['N_PROF']
        z_ipar= np.zeros(N_PROF) # Quenching layer depth

        for i in range(N_PROF):
            if ~np.isnan(np.nansum(ds.CHLA_ADJUSTED[i,:])):
                x=ds.LIGHT_ADJUSTED[i,:].values - iPAR
                x = x[~np.isnan(x)]
                idx_ipar= np.abs((x)).argmin(axis=0)
                z_ipar[i]=ds.N_LEVELS[idx_ipar].values
                if (z_ipar[i]< ds.MLD[i]) :
                    for j in range(len(ds.CHLA_ADJUSTED[i,:idx_ipar].values)):
                        ds['CHLA_CORRECTED'][i,j] = ds.CHLA_ADJUSTED[i,idx_ipar]

        ds['ipar15'] = xr.DataArray(z_ipar,dims=['N_PROF'])