#
# Copyright 2018-2020 École Polytechnique Fédérale de Lausanne (EPFL) and
# Paul Scherrer Institut (PSI).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import pandas as pd
import numpy as np
import datetime
from scipy.interpolate import interp1d # for afc correction
from pathlib import Path

from pyantarctica import aceairsea


def ensure_tz_UTC(df):
    """
        Function to ensure that the tzinfo of a data frame is in UTC

        :param df: a pandas data frame with datetime_index
        :returns: the dataframe with df.index.tzinfo='UTC'
    """
    # ensure tzinfo is in UTC
    if (str(df.index.tzinfo)=='None'):
        #print('The time stamp is provided is tz naive. It is assumed that the times provided are in UTC')
        df=df.tz_localize(tz='UTC')
    else:
        #print('The provided time zone of the data is '+str(df.index.tzinfo)+'. This is now converted to a tz aware time stampe in UTC')
        df.index=df.index.tz_convert('UTC')
    return df

def read_assimilation_list():
    """
        Function to read the list of time stamps where data from UBXH3 was asssimilated

        :returns: a dataframe with the assimilation time stamps
    """
    # TODO change path
    UBXH3_file = Path('./data/ubxh3_assimilation_list/77040_UBXH3_2017012618_2017021818.txt')
    # read UTC, latitude, longitude of occations where station UBXH3 was assimilated into the IFS
    #UBXH3_assimilated = pd.read_csv(UBXH3_file, header=None, delimiter=' +', parse_dates={'date_time':[1, 2]}, names=['station',  'yyyymmdd', 'hh', 'latitude', 'longitude'], engine='python')
    UBXH3_assimilated = pd.read_csv(UBXH3_file, delimiter=' +', parse_dates={'date_time':[1, 2]}, engine='python')
    UBXH3_assimilated.set_index( (pd.to_datetime(UBXH3_assimilated.date_time, format="%Y-%m-%d %H:%M:%S")), drop=True, inplace=True )
    UBXH3_assimilated.drop(columns=['date_time'], inplace=True)
    UBXH3_assimilated=UBXH3_assimilated.tz_localize(tz='UTC')
    return UBXH3_assimilated



def read_distance2land():
    """
        Function to read the ship's distance to land

        :returns: a dataframe with the distance to land information
    """
    # TODO change path
    # NOTE UNITS OF distance are in meter!
    dist2land_file = Path('../cruise-track-distance-to-land/data/distance_to_land_corrected_on_small_islands/dist_to_land_incl_small_islands.csv')
    dist2land_file = Path('./data/distance_to_the_nearest_10/dist_to_land_incl_small_islands.csv')
    dist2land = pd.read_csv(dist2land_file)
    if ('timest_' in dist2land.columns):
        dist2land.rename(columns={'timest_': 'date_time'}, inplace=True)
    dist2land.set_index(pd.DatetimeIndex(dist2land.date_time), inplace=True) 
    dist2land.drop(columns=['date_time'], inplace=True)
    
    dist2land = ensure_tz_UTC(dist2land)
    
    return dist2land

def read_era5_data():
    """
        Function to read the ERA-5 reanalysis data which was interpolated onto the ship's track.
        Calculates the Obukov length scale

        :returns: a dataframe with the interpolated ERA-5 reanalysis data
    """
    era5_csv_file = Path('../ecmwf-interpolation-to-cruise-track/data/ecmwf-era5-interpolated-to-cruise-track/era5-on-cruise-track-5min-legs0-4-nearest.csv')
    era5_csv_file = Path('./data/era5_reanalysis_results_10/era5-on-cruise-track-5min-legs0-4-nearest.csv')
    era5 = pd.read_csv(era5_csv_file)
    if 1:
        era5.rename(columns={
            '10m_u_component_of_neutral_wind':'u10n',
            '10m_v_component_of_neutral_wind':'v10n',
            '10m_u_component_of_wind':'u10',
            '10m_v_component_of_wind':'v10',
            'air_density_over_the_oceans':'p140209',
            '2m_dewpoint_temperature':'d2m',
            '2m_temperature':'t2m',
            'boundary_layer_height':'blh',
            'cloud_base_height':'cbh',
            'convective_precipitation':'cp',
            'convective_snowfall':'csf',
            'friction_velocity':'zust',
            'high_cloud_cover':'hcc',
            'land_sea_mask':'lsm',
            'large_scale_precipitation':'lsp',
            'large_scale_snowfall':'lsf',
            'low_cloud_cover':'lcc',
            #'mdts':'mean_direction_of_total_swell',
            #'mdww':'mean_direction_of_wind_waves',
            #'msqs':'mean_square_slope_of_waves',
            'mean_surface_latent_heat_flux':'mslhf',
            'mean_surface_sensible_heat_flux':'msshf',
            #'mcc':'medium_cloud_cover',
            #'pp1d':'peak_wave_period',
            #'tmax':'period_corresponding_to_maximum_individual_wave_height',
            'sea_ice_cover':'siconc',
            'sea_surface_temperature':'sst',
            #'shts':'significant_height_of_total_swell',
            #'shww':'significant_height_of_wind_waves',
            'skin_temperature':'skt',
            #'tcc':'total_cloud_cover',
            #'tciw':'total_column_cloud_ice_water',
            #'tclw':'total_column_cloud_liquid_water'
        }, inplace=True)

    if ('timest_' in era5.columns):
        era5.rename(columns={'timest_': 'date_time'}, inplace=True)
    era5.set_index(pd.DatetimeIndex(era5.date_time), inplace=True) 
    era5.drop(columns=['date_time'], inplace=True)

    era5 = ensure_tz_UTC(era5)
    
    era5['WS10']=np.sqrt( np.square(era5.u10)+np.square(era5.v10) ) # 10 meter wind speed (U_10)
    era5['WS10N']=np.sqrt( np.square(era5.u10n)+np.square(era5.v10n) ) # 10 meter NEUTRAL wind speed (U_10N)

    era5 = era5.assign(msshf_original = era5['msshf']) # copy of the orignal sensible heaf flux
    # here we fix the values of msshf near land to follow the bulk sea water behaviour
    era5['sst'][0]=era5['sst'][np.where(era5['sst']>0)[0][0]] # first 8 entries are NaN need to fill first value for interpolation to work
    era5['sst'] = era5['sst'].interpolate() # fill missing values with basic interpolation. This is good enough
    era5['msshf'][era5['lsm']>0] = (1.7*(era5['t2m']-era5['sst'])*era5['WS10N'])[era5['lsm']>0]

    era5['LMO'] = aceairsea.LMoninObukov_bulk(era5['WS10N'],-era5.msshf,-era5.mslhf,(era5.t2m-273.15)) # Monin-Obukov Length Scale
    era5['ustar'] = aceairsea.coare_u2ustar (era5['WS10N'], input_string='u2ustar', coare_version='coare3.5', TairC=(era5.t2m-273.15), z=10, zeta=0)

    era5['LSM']=era5['lsm']
    era5['SIF']=era5['siconc']
    era5.at[np.isnan(era5['SIF']),'SIF']=0

    #'T2M', 'LMO', 'SKT',
    era5['T2M']=era5['t2m']-273.15 # air temperatue in Celsius
    
    return era5



def read_and_filter_wind_data():
    """
        Function to read the raw wind data from all 5 legs

        :returns: a dataframe with the combined wind data with spurious observations removed
    """
    wind_csv_file_folder = './data/summary_raw_wind_data_fr_12/'
    wind_csv_file_name0 = 'metdata_wind_20161117_20161216.csv'
    wind_csv_file_name1 = 'metdata_wind_20161220_20170118.csv'
    wind_csv_file_name2 = 'metdata_wind_20170122_20170223.csv'
    wind_csv_file_name3 = 'metdata_wind_20170226_20170319.csv'
    wind_csv_file_name4 = 'metdata_wind_20170322_20170411.csv'

    df_wind0 = pd.read_csv(wind_csv_file_folder+wind_csv_file_name0)
    df_wind1 = pd.read_csv(wind_csv_file_folder+wind_csv_file_name1)
    df_wind2 = pd.read_csv(wind_csv_file_folder+wind_csv_file_name2)
    df_wind3 = pd.read_csv(wind_csv_file_folder+wind_csv_file_name3)
    df_wind4 = pd.read_csv(wind_csv_file_folder+wind_csv_file_name4)
    
    df_wind0 = df_wind0.set_index( (pd.to_datetime(df_wind0.date_time, format="%Y-%m-%d %H:%M:%S")+pd.to_timedelta(df_wind0.TIMEDIFF, unit='s'))  ) # assing the time stamp
    
    df_wind1 = df_wind1.set_index( (pd.to_datetime(df_wind1.date_time, format="%Y-%m-%d %H:%M:%S")+pd.to_timedelta(df_wind1.TIMEDIFF, unit='s'))  ) # assing the time stamp
    
    df_wind2 = df_wind2.set_index( (pd.to_datetime(df_wind2.date_time, format="%Y-%m-%d %H:%M:%S")+pd.to_timedelta(df_wind2.TIMEDIFF, unit='s'))  ) # assing the time stamp
    
    df_wind3 = df_wind3.set_index( (pd.to_datetime(df_wind3.date_time, format="%Y-%m-%d %H:%M:%S")+pd.to_timedelta(df_wind3.TIMEDIFF, unit='s'))  ) # assing the time stamp
    
    df_wind4 = df_wind4.set_index( (pd.to_datetime(df_wind4.date_time, format="%Y-%m-%d %H:%M:%S")+pd.to_timedelta(df_wind4.TIMEDIFF, unit='s'))  ) # assing the time stamp
    
    frames = [df_wind0, df_wind1, df_wind2, df_wind3, df_wind4] # concatenate the 4 wind data files
    df_wind = pd.concat(frames)
       
    max_wind_WS = 40    # True wind speed maximum value (arbitrary)
    max_wind_WSR = 100 # bad data flaged with 999 # one spike in WSR2 at 700 
    # removes all bad directions (out of 0-360) for 1 but not for 2!

    print('WSR1: ' +str(np.sum(df_wind.WSR1>-1))+' samples for WSR1')
    print('WSR2: ' +str(np.sum(df_wind.WSR2>-1))+' samples for WSR2')

    
    print('Removing ' +str(np.sum(df_wind.WSR1>max_wind_WSR))+' samples for WSR1 > 100 m/s')

    print('Removing ' +str(np.sum(df_wind.WSR2>max_wind_WSR))+' samples for WSR2 > 100 m/s')

    df_wind.at[df_wind.WSR1>max_wind_WSR, 'WD1'] = np.nan
    df_wind.at[df_wind.WSR1>max_wind_WSR, 'WS1'] = np.nan
    df_wind.at[df_wind.WSR1>max_wind_WSR, 'WDR1'] = np.nan
    df_wind.at[df_wind.WSR1>max_wind_WSR, 'WSR1'] = np.nan

    df_wind.at[df_wind.WS1>max_wind_WS, 'WD1'] = np.nan
    df_wind.at[df_wind.WS1>max_wind_WS, 'WS1'] = np.nan

    df_wind.at[df_wind.WSR2>max_wind_WSR, 'WD2'] = np.nan
    df_wind.at[df_wind.WSR2>max_wind_WSR, 'WS2'] = np.nan
    df_wind.at[df_wind.WSR2>max_wind_WSR, 'WDR2'] = np.nan
    df_wind.at[df_wind.WSR2>max_wind_WSR, 'WSR2'] = np.nan

    df_wind.at[df_wind.WS2>max_wind_WS, 'WD2'] = np.nan
    df_wind.at[df_wind.WS2>max_wind_WS, 'WS2'] = np.nan
    
    # remove all W2 readings for WDR2>360
    print('Removing ' +str(np.sum(df_wind.WDR1>360))+' samples for WDR1 > 360deg')
    print('Removing ' +str(np.sum(df_wind.WDR2>360))+' samples for WDR2 > 360deg')

    df_wind.at[df_wind.WDR2>360, 'WD2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WS2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WSR2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WD2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WDR2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WSR2'] = np.nan
    df_wind.at[df_wind.WDR2>360, 'WS2'] = np.nan    
    df_wind.at[df_wind.WDR2>360, 'WDR2'] = np.nan

    dWSR = np.abs(df_wind.WSR1-df_wind.WSR2) # difference between the relative wind speeds
    print('Removing ' +str(np.sum(dWSR>20))+' samples for dWSR > 20 m/s')
    
    df_wind.at[dWSR>20, 'WD1'] = np.nan
    df_wind.at[dWSR>20, 'WS1'] = np.nan
    df_wind.at[dWSR>20, 'WDR1'] = np.nan
    df_wind.at[dWSR>20, 'WSR1'] = np.nan
    df_wind.at[dWSR>20, 'WD2'] = np.nan
    df_wind.at[dWSR>20, 'WS2'] = np.nan
    df_wind.at[dWSR>20, 'WDR2'] = np.nan
    df_wind.at[dWSR>20, 'WSR2'] = np.nan

    # S2 hase some biased low wind speed for 6 hours -> remove WSR2 and WDR2 for these
    print('removing samples for s2 2017-04-05 06:00:00 till 2017-04-05 12:00:00')
    df_wind.at[ ((df_wind.index>pd.to_datetime("2017-04-05 06:00:00", format="%Y-%m-%d %H:%M:%S").tz_localize(tz='UTC'))&(df_wind.index<pd.to_datetime("2017-04-05 12:00:00", format="%Y-%m-%d %H:%M:%S").tz_localize(tz='UTC')))==1, 'WSR2'] = np.NaN
    df_wind.at[ ((df_wind.index>pd.to_datetime("2017-04-05 06:00:00", format="%Y-%m-%d %H:%M:%S").tz_localize(tz='UTC'))&(df_wind.index<pd.to_datetime("2017-04-05 12:00:00", format="%Y-%m-%d %H:%M:%S").tz_localize(tz='UTC')))==1, 'WDR2'] = np.NaN
    
    # remove extraneous columns
    df_wind = df_wind.drop(columns=['COG', 'CLOUDTEXT', 'date_time'])
    
    df_wind.index.name = 'date_time'
    
    return df_wind

def read_minute_ship_track():
    """
        Function to read the 1-minute resolution ship track velocity file

        :returns: a dataframe with the data
    """
    gps_csv_file_folder = './data/oneminute_average_cruis_11/'
    gps_csv_file_name = 'cruise-track-1min-legs0-4.csv'
    df_gps = pd.read_csv(gps_csv_file_folder+gps_csv_file_name)
    df_gps = df_gps.set_index(pd.to_datetime(df_gps.date_time, format="%Y-%m-%dT%H:%M:%S")) # assing the time stamp
    df_gps.drop(columns=['date_time'], inplace=True)

    df_gps = df_gps.rename(columns={'platform_speed_wrt_sea_water_east':'velEast', 
                                      'platform_speed_wrt_sea_water_north':'velNorth', 
                                      'platform_orientation':'HEADING', 
                                      'platform_course':'COG', 
                                      'platform_speed_wrt_ground':'SOG'})
    
    return df_gps

def afc_expand(df):
    """
        Function to calculate U10N from wind speed measured at height Zanemometer and accounting for atmospheric stability

        :param df: data fram needs to contain a variable df.R denoting direction over the interval [-180 180)
        :returns: df: data fram expanded over the variable df.R
    """
    df_left = df.copy()
    df_left.R=df_left.R-360
    df_right = df.copy()
    df_right.R=df_right.R+360
    df = pd.concat([df_left, df, df_right])
    return df

def wind_merge_gps_afc_option(afc_correction_factor_files=[]):
    """
        Function to merge the wind speed and ship track velocity data at 1 minute resolution
        Optional the air flow distortion correction of the observed relative wind speed and direction is performed prior to the motion correction

        :returns: a dataframe with the combined wind (true and relative) and ship velocity data
    """
    merge_at_nseconds = 60
    df_gps = read_minute_ship_track();  print("read_minute_ship_track")

    df_wind = read_and_filter_wind_data();  print("read_and_filter_wind_data ... ")
    
    if len(afc_correction_factor_files)>0: # if afc data provided
        # interpolate onto wind series and calculate afc wind speed and direction

        afc_loess1 = pd.read_csv(afc_correction_factor_files+'1.csv', index_col=False)
        afc_loess2 = pd.read_csv(afc_correction_factor_files+'2.csv', index_col=False)
        if 'Unnamed: 0' in afc_loess1.columns:
            afc_loess1=afc_loess1.drop(columns={'Unnamed: 0'})
        if 'Unnamed: 0' in afc_loess2.columns:
            afc_loess2=afc_loess2.drop(columns={'Unnamed: 0'})
        
        if "R" not in afc_loess1.columns:
            afc_loess1.rename(columns={
                'wind_from_direction_relative_to_platform': 'R',
                'mean_of_wind_speed_bias': 'A',
                'mean_of_wind_direction_bias': 'D',
                'median_of_wind_speed_bias': 'A_median',
                'median_of_wind_direction_bias': 'D_median',
                'uncertainty_of_wind_speed_bias': 'A_err',
                'uncertainty_of_wind_direction_bias': 'D_err',
                'number_of_samples': 'samples'
            }, inplace=True)
            
        if "R" not in afc_loess2.columns:
            afc_loess2.rename(columns={
                'wind_from_direction_relative_to_platform': 'R',
                'mean_of_wind_speed_bias': 'A',
                'mean_of_wind_direction_bias': 'D',
                'median_of_wind_speed_bias': 'A_median',
                'median_of_wind_direction_bias': 'D_median',
                'uncertainty_of_wind_speed_bias': 'A_err',
                'uncertainty_of_wind_direction_bias': 'D_err',
                'number_of_samples': 'samples'
            }, inplace=True)
        
        afc_loess1=afc_expand(afc_loess1)
        afc_loess2=afc_expand(afc_loess2)

        # interpolate bias onto time series
        r1=((df_wind.WDR1)+180)%360-180
        fA1 = interp1d(afc_loess1.R, afc_loess1.A, kind='linear')
        a1 = fA1(r1)
        fD1 = interp1d(afc_loess1.R, afc_loess1.D, kind='linear')
        d1 = fD1(r1)
        # interpolate bias onto time series
        r2=((df_wind.WDR2)+180)%360-180
        fA2 = interp1d(afc_loess2.R, afc_loess2.A, kind='linear')
        a2 = fA2(r2)
        fD2 = interp1d(afc_loess2.R, afc_loess2.D, kind='linear')
        d2 = fD2(r2)

        # WSR = WSR/a
        # WDR =  WDR-d
        #s1
        df_wind.WSR1 = np.true_divide(df_wind.WSR1,a1)
        df_wind.WDR1 = (df_wind.WDR1-d1)
        #s2
        df_wind.WSR2 = np.true_divide(df_wind.WSR2,a2)
        df_wind.WDR2 = (df_wind.WDR2-d2)

    # calculate uR,vR and u,v from SR, DR and WS,WR
    # uR are positive for wind along ships main axis
    # vR are positive for wind in port direction
    # u are positive for wind blowing East
    # v are positive for wind blowing North
    df_wind = df_wind.assign(uR1=df_wind.WSR1*np.cos(np.deg2rad(180-df_wind.WDR1)))
    df_wind = df_wind.assign(vR1=df_wind.WSR1*np.sin(np.deg2rad(180-df_wind.WDR1)))
    df_wind = df_wind.assign(u1=df_wind.WS1*np.cos(np.deg2rad(270-df_wind.WD1)))
    df_wind = df_wind.assign(v1=df_wind.WS1*np.sin(np.deg2rad(270-df_wind.WD1)))
    #
    df_wind = df_wind.assign(uR2=df_wind.WSR2*np.cos(np.deg2rad(180-df_wind.WDR2)))
    df_wind = df_wind.assign(vR2=df_wind.WSR2*np.sin(np.deg2rad(180-df_wind.WDR2)))
    df_wind = df_wind.assign(u2=df_wind.WS2*np.cos(np.deg2rad(270-df_wind.WD2)))
    df_wind = df_wind.assign(v2=df_wind.WS2*np.sin(np.deg2rad(270-df_wind.WD2)))

    # add heading sin cos for proper averaging
    df_wind = df_wind.assign(hdg_cos=np.cos(np.deg2rad(df_wind.HEADING)))
    df_wind = df_wind.assign(hdg_sin=np.sin(np.deg2rad(df_wind.HEADING)))

    # assing apparent wind in eart,north coordinates
    df_wind = df_wind.assign(uA1=df_wind.WSR1*np.cos(np.deg2rad(270-df_wind.HEADING-df_wind.WDR1)))
    df_wind = df_wind.assign(vA1=df_wind.WSR1*np.sin(np.deg2rad(270-df_wind.HEADING-df_wind.WDR1)))
    df_wind = df_wind.assign(uA2=df_wind.WSR2*np.cos(np.deg2rad(270-df_wind.HEADING-df_wind.WDR2)))
    df_wind = df_wind.assign(vA2=df_wind.WSR2*np.sin(np.deg2rad(270-df_wind.HEADING-df_wind.WDR2)))



    nSample_min = 6 # out of 20 possible (usually 12+/2)

    nWSR1=df_wind.WSR1.resample(str(merge_at_nseconds)+'S', loffset = datetime.timedelta(seconds=(merge_at_nseconds*0.5))).count() # calculate 6 second average
    nWSR2=df_wind.WSR2.resample(str(merge_at_nseconds)+'S', loffset = datetime.timedelta(seconds=(merge_at_nseconds*0.5))).count() # calculate 6 second average
    #plt.plot(nWSR1,'.')
    # strangely some bins contain more than nSec/3 wind samples!!!

    print( 'WSR1 removing ' + str( np.sum((nWSR1>0) & (nWSR1<nSample_min)) ) + ' of '+str(np.sum(nWSR1>0)) )
    print( 'WSR2 removing ' + str( np.sum((nWSR2>0) & (nWSR2<nSample_min)) ) + ' of '+str(np.sum(nWSR1>0)) )

    if 1: # compute HEADING DIFF accross 1min interval
        HEADING = df_wind.HEADING.copy();
        HEADING_MAX=HEADING.resample(str(merge_at_nseconds)+'S').max()
        HEADING_MIN=HEADING.resample(str(merge_at_nseconds)+'S').min()
        HEADING = (HEADING-180)%360
        HEADING_MAX_=HEADING.resample(str(merge_at_nseconds)+'S').max()
        HEADING_MIN_=HEADING.resample(str(merge_at_nseconds)+'S').min()
        HEADING_DIFF = np.min([(HEADING_MAX-HEADING_MIN), (HEADING_MAX_-HEADING_MIN_)], axis=0)



    df_wind=df_wind.resample(str(merge_at_nseconds)+'S', loffset = datetime.timedelta(seconds=(merge_at_nseconds*0.5))).mean() # calculate 6 second average

    if 1:
        df_wind = df_wind.assign(HEADING_DIFF=HEADING_DIFF)

        df_wind = df_wind.assign(nWSR1=nWSR1)
        df_wind = df_wind.assign(nWSR2=nWSR2)

    # set data with too few readings to NaN
    df_wind.at[nWSR1<nSample_min, 'uR1'] = np.nan
    df_wind.at[nWSR1<nSample_min, 'vR1'] = np.nan
    df_wind.at[nWSR1<nSample_min, 'uA1'] = np.nan
    df_wind.at[nWSR1<nSample_min, 'vA1'] = np.nan
    df_wind.at[nWSR2<nSample_min, 'uR2'] = np.nan
    df_wind.at[nWSR2<nSample_min, 'vR2'] = np.nan
    df_wind.at[nWSR2<nSample_min, 'uA2'] = np.nan
    df_wind.at[nWSR2<nSample_min, 'vA2'] = np.nan

    # rebuild the angles:
    df_wind.HEADING = np.rad2deg(np.arctan2(df_wind.hdg_sin, df_wind.hdg_cos)) % 360
    df_wind.WD1 = (270 - np.rad2deg(np.arctan2(df_wind.v1, df_wind.u1)) )% 360
    df_wind.WDR1 = (180 - np.rad2deg(np.arctan2(df_wind.vR1, df_wind.uR1)) )% 360
    df_wind.WD2 = (270 - np.rad2deg(np.arctan2(df_wind.v2, df_wind.u2)) )% 360
    df_wind.WDR2 = (180 - np.rad2deg(np.arctan2(df_wind.vR2, df_wind.uR2)) )% 360
    # recalcualte the speeds as vector average
    df_wind.WS1 = np.sqrt( np.square(df_wind.v1) + np.square(df_wind.u1) )
    df_wind.WS2 = np.sqrt( np.square(df_wind.v2) + np.square(df_wind.u2) )
    df_wind.WSR1 = np.sqrt( np.square(df_wind.vR1) + np.square(df_wind.uR1) )
    df_wind.WSR2 = np.sqrt( np.square(df_wind.vR2) + np.square(df_wind.uR2) )

    
    if ('HEADING' in df_gps.columns):
        df_gps = df_gps.drop(columns=['HEADING']) # avoid to have the same fields in both data frames
    
    # merge with GPS on 60sec basis
    df_wind = df_wind.merge(df_gps, left_on='date_time', right_on='date_time', how='left')

    df_wind.COG = (90-np.rad2deg(np.arctan2(df_wind.velNorth,df_wind.velEast))) % 360 # recompute COG from averaged North/Easte velocities


    # True wind correction using GPS data

    # use average apparent wind for motion correction
    df_wind.u1 = df_wind.uA1 + df_wind.velEast
    df_wind.v1 = df_wind.vA1 + df_wind.velNorth
    df_wind.u2 = df_wind.uA2 + df_wind.velEast
    df_wind.v2 = df_wind.vA2 + df_wind.velNorth

    # rebuild the angles:
    df_wind.WD1 = (270 - np.rad2deg(np.arctan2(df_wind.v1, df_wind.u1)) )% 360
    df_wind.WD2 = (270 - np.rad2deg(np.arctan2(df_wind.v2, df_wind.u2)) )% 360
    # recalcualte the speeds as vector average
    df_wind.WS1 = np.sqrt( np.square(df_wind.v1) + np.square(df_wind.u1) )
    df_wind.WS2 = np.sqrt( np.square(df_wind.v2) + np.square(df_wind.u2) )

    return df_wind