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
from sklearn.utils import resample # for bootstrapping


from collections import defaultdict

from pathlib import Path

from pyantarctica import windvectorcoordinates as wvc
from pyantarctica import aceairsea as aceairsea 
from pyantarctica import dataset as dataset

from pyantarctica.datafilter import outliers_iqr_noise


def Rdir_bin_intervals(R1,A1,QC,binl,binu,min_in_bin=12,find_IQR_outlier=True,NOISE=0.1, Weights_a=[], BOOTSTRAP=False):
    """
        Function to bin data over specified direction bins
        with Options to 
        - identify and exclude outliers
        - calculate weighted averages
        - apply bootstrapping to estimate the uncertainty of the mean

        :param R1: independent variable (e.g. relative wind direction)
        :param A1: dependent variable (e.g. wind-speed-ratio)
        :param QC: quality flag (only data A1[QC==True] are used)
        :param binl: list of the lower limits of the bins
        :param binu: list of the upper limits of the bins
        :param min_in_bin: Minimum number of observations per bin required to return a non-NaN result
        :param find_IQR_outlier: True/False, if True, data for each bin are screened for outlieres using outliers_iqr_noise(A1,NOISE)
        :param NOISE: noise level of A1 used in outliers_iqr_noise(A1,NOISE) to avoid rejection of digital noise
        :param Weights_a: (OPTIONAL) Weights of A1 to be used in the averaging
        :param BOOTSTRAP: True/Flase, if True the uncertainty is estimated via bootstrapping, if False, normal standard deviations are calcualted
        :returns: bin_data: DICT containing the bin averages and uncertainties of R1 and A1
        :returns: outlier_mask: True/False list (length of R1) identifying the detected outliers
        :returns: outlier_count:
    """
    # Weights used to better estimate the mean
    if len(Weights_a)==0:
        Weights_a = np.ones_like(R1)
    Weights_a[np.isnan(Weights_a)]=0 # set indefinite weights to zero
    
    bin_data = defaultdict(list)
    outlier_mask = np.ones_like(QC) #
    outlier_count = np.zeros_like(QC) # count how often a point is flagged as outlier
    binl=wvc.ang180(binl)
    binu=wvc.ang180(binu)
    R1=wvc.ang180(R1)
    
    for jbin in np.arange(0,len(binl),1):
        if binl[jbin]<binu[jbin]:
            in_bin = (((R1>=binl[jbin]) & (R1<binu[jbin])) & QC )
        elif binl[jbin]>binu[jbin]: # if the bin goes accross 180
            in_bin = (((R1>=binl[jbin]) | (R1<binu[jbin])) & QC )
        else:
            print('bin edges must not be equal!'); return;
        
        if sum(in_bin)>min_in_bin:
            if find_IQR_outlier:
                outliers = outliers_iqr_noise(A1[in_bin==1],NOISE) # minimu deviation of 5%  or 10% ?
            else:
                outliers=[]
            outlier_mask[np.where(in_bin)[0][outliers]]=0
            outlier_count[np.where(in_bin)[0][outliers]]=(outlier_count[np.where(in_bin)[0][outliers]]+1)

            if binl[jbin]<binu[jbin]:
                in_bin = (((R1>=binl[jbin]) & (R1<binu[jbin])) & QC & outlier_mask)
            elif binl[jbin]>binu[jbin]: # if the bin goes accross 180
                in_bin = (((R1>=binl[jbin]) | (R1<binu[jbin])) & QC & outlier_mask)
            
            if sum(in_bin)>min_in_bin:
                binStd = np.nanstd(A1[in_bin==1])
                binMedian = np.nanmedian(A1[in_bin==1])
                binMean = np.nanmean(A1[in_bin==1])
                binMean = np.average(A1[in_bin==1], weights=Weights_a[in_bin==1])
                binSum = sum(in_bin)
                
                

                if len(Weights_a)==np.sum(Weights_a):
                    binErr = binStd/np.sqrt(binSum)
                else:
                    binErr = 1/np.sqrt(np.sum(Weights_a[in_bin==1]))
                    # https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
                    
                    
                if BOOTSTRAP:
                    N_iter = 100
                    binMeans = np.arange(0,N_iter)*np.NaN
                    for j in np.arange(0,N_iter):
                        A1_rs, weights_rs = resample( A1[in_bin==1],Weights_a[in_bin==1], random_state=j+42)
                        binMeans[j]=np.average(A1_rs, weights=weights_rs)
                        binErr = np.std(binMeans)
                
                    
            else:
                binMean=np.nan; binMedian=np.nan; binStd=np.nan; binSum=sum(in_bin); binErr=np.nan; binErr2=np.nan
        else:
            binMean=np.nan; binMedian=np.nan; binStd=np.nan; binSum=sum(in_bin); binErr=np.nan; binErr2=np.nan
            
        
        if binl[jbin]<binu[jbin]:
            bin_data['x_mean'].append(np.nanmean(R1[in_bin==1]))
        elif binl[jbin]>binu[jbin]: 
            # if the bin goes accross 180
            # need to bring lower part on top side, average and then ang180 again
            r1_Leq180 = R1[(((R1>=binl[jbin]) ) & QC & outlier_mask)]
            r1_Geq180 = R1[(((R1<binu[jbin]) ) & QC & outlier_mask)]
            bin_data['x_mean'].append(wvc.ang180  ( np.nanmean( np.concatenate([r1_Leq180, (r1_Geq180+360) ]) ) ) )
            
        bin_data['y_mean'].append(binMean)
        bin_data['y_median'].append(binMedian)
        bin_data['y_std'].append(binStd)
        bin_data['samples'].append(binSum)
        bin_data['y_err'].append(binErr)


    bin_data = pd.DataFrame(bin_data)
    bin_data['x_mean']=wvc.ang180(bin_data['x_mean'])
    bin_data = bin_data.sort_values('x_mean') # sort along R1 direction
    return bin_data, outlier_mask, outlier_count

def tryoshnikov_afc(t,d,s,D,S, QC, high_res=False,find_IQR_outlier=False, BOOTSTRAP=False, Weights_a=[], Weights_d=[]):
    """
        Function to run the flow distortion anlaysis for the ACE data
        with Options to 
        - identify and exclude outliers
        - calculate weighted averages (by specifying the weights of the observations)
        - apply bootstrapping to estimate the uncertainty of the mean
        - use disjunct or overlapping wind direction intervals

        :param t: date time index
        :param d: observed relative wind direction
        :param s: observed relative wind speed
        :param D: expected relative wind direction (calculated from the freestreem reference)
        :param S: expected relative wind speed (calculated from the freestreem reference)
        :param QC: quality flag (only data [QC==True] are used)
        :high_res: True/False: if True, the wind direction bins are overlapping, if false, the wind direction bins are disjunct
        :param find_IQR_outlier: True/False, if True, data for each bin are screened for outlieres using outliers_iqr_noise(A1,NOISE)
        :param NOISE: noise level of A1 used in outliers_iqr_noise(A1,NOISE) to avoid rejection of digital noise
        :param BOOTSTRAP: True/Flase, if True the uncertainty is estimated via bootstrapping, if False, normal standard deviations are calcualted
        :param Weights_a: (OPTIONAL) Weights for s/S to be used in the averaging
        :param Weights_d: (OPTIONAL) Weights for d-D to be used in the averaging
        :returns: radqc_: Dataframe containing R=ang180(d),A=s/S,D=d-D,QC=QC,outliers={outliers identified via IQR filter} 
        :returns: afc_: Dataframe containing the bin mean, median and std values of R, A, D,
    """
    # run ratio and dirdiff analysis
    # optional filter based on IQR
    # winddirection bins are hardset, adapted to the dataset
    min_in_bin = 7

    R = wvc.ang180(d)
    A = s/S
    dD = wvc.ang180(wvc.ang180(d)-wvc.ang180(D))
    if high_res:
        binl = np.concatenate( [np.arange(-180,-130,2), np.arange(-130,180,1)] ) # was at -140
        binu = np.concatenate( [np.arange(-180,-130,2)+10, np.arange(-130,180,1)+5] )
    else:
        binl = np.concatenate( [np.arange(-180,-130,10), np.arange(-130,185,5)] ) # was at -140
        binu = np.concatenate( [np.arange(-180,-130,10)+10, np.arange(-130,185,5)+5] )

    outliers = np.zeros_like(QC)

    bin_data_A, outlier_mask_A, outlier_count_A = Rdir_bin_intervals(R,A,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.05, Weights_a=Weights_a, BOOTSTRAP=False)
    bin_data_D, outlier_mask_D, outlier_count_D = Rdir_bin_intervals(R,dD,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.1, Weights_a=Weights_d, BOOTSTRAP=False)

    outliers = outliers | (outlier_count_A | outlier_count_A) # collect all outliers
    if find_IQR_outlier:
        QC = (QC & outlier_mask_A & outlier_mask_D)
    if 1: # always run a second time
        bin_data_A, outlier_mask_A, outlier_count_A = Rdir_bin_intervals(R,A,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.05, Weights_a=Weights_a, BOOTSTRAP=BOOTSTRAP)
        bin_data_D, outlier_mask_D, outlier_count_D = Rdir_bin_intervals(R,dD,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.1, Weights_a=Weights_d, BOOTSTRAP=BOOTSTRAP)
        QC = QC & (outlier_mask_A & outlier_mask_D)
        outliers = outliers | (outlier_count_A | outlier_count_A) # collect all outliers


    afc_ = pd.DataFrame({'R':bin_data_A['x_mean'],'A':bin_data_A['y_mean'],'D':bin_data_D['y_mean'], 
                         'A_median':bin_data_A['y_median'],'D_median':bin_data_D['y_median'],
                         'A_err':(bin_data_A['y_err']), 'D_err':(bin_data_D['y_err']), 'samples':bin_data_A['samples']})
    # afc data
    radqc_ = pd.DataFrame({'R':R,'A':A,'D':dD,'QC':QC,'outliers':outliers})
    

    
    return radqc_, afc_


def Rdir_bin_intervals_unique(R1,A1,S0,QC,binl,binu,min_in_bin=12,find_IQR_outlier=True,NOISE=0.1, Weights_a=[], BOOTSTRAP=False):
    """
        Function to bin data over specified direction bins, identifying if the reference observations are not independent.
        with Options to 
        - identify and exclude outliers
        - calculate weighted averages
        - apply bootstrapping to estimate the uncertainty of the mean

        :param R1: independent variable (e.g. relative wind direction)
        :param A1: dependent variable (e.g. wind-speed-ratio)
        :param S0: variable used to identify if the observations are truely unique per wind direction sector (e.g. use ERA-5 U10N)
        :param QC: quality flag (only data A1[QC==True] are used)
        :param binl: list of the lower limits of the bins
        :param binu: list of the upper limits of the bins
        :param min_in_bin: Minimum number of observations per bin required to return a non-NaN result
        :param find_IQR_outlier: True/False, if True, data for each bin are screened for outlieres using outliers_iqr_noise(A1,NOISE)
        :param NOISE: noise level of A1 used in outliers_iqr_noise(A1,NOISE) to avoid rejection of digital noise
        :param Weights_a: (OPTIONAL) Weights of A1 to be used in the averaging
        :param BOOTSTRAP: True/Flase, if True the uncertainty is estimated via bootstrapping, if False, normal standard deviations are calcualted
        :returns: bin_data: DICT containing the bin averages and uncertainties of R1 and A1
        :returns: outlier_mask: True/False list (length of R1) identifying the detected outliers
        :returns: outlier_count:
    """ 
    # Weights used to better estimate the mean
    if len(Weights_a)==0:
        Weights_a = np.ones_like(R1)
    Weights_a[np.isnan(Weights_a)]=0 # set indefinite weights to zero
    
    bin_data = defaultdict(list)
    outlier_mask = np.ones_like(QC) #
    outlier_count = np.zeros_like(QC) # count how often a point is flagged as outlier
    binl=wvc.ang180(binl)
    binu=wvc.ang180(binu)
    R1=wvc.ang180(R1)
    
    for jbin in np.arange(0,len(binl),1):
        if binl[jbin]<binu[jbin]:
            in_bin = (((R1>=binl[jbin]) & (R1<binu[jbin])) & QC )
        elif binl[jbin]>binu[jbin]: # if the bin goes accross 180
            in_bin = (((R1>=binl[jbin]) | (R1<binu[jbin])) & QC )
        else:
            print('bin edges must not be equal!'); return;
        
        #if sum(in_bin)>min_in_bin:
        if len(np.unique(S0[in_bin]))>min_in_bin:
            if find_IQR_outlier:
                outliers = outliers_iqr_noise(A1[in_bin==1],NOISE) # minimu deviation of 5%  or 10% ?
            else:
                outliers=[]
            outlier_mask[np.where(in_bin)[0][outliers]]=0
            outlier_count[np.where(in_bin)[0][outliers]]=(outlier_count[np.where(in_bin)[0][outliers]]+1)

            if binl[jbin]<binu[jbin]:
                in_bin = (((R1>=binl[jbin]) & (R1<binu[jbin])) & QC & outlier_mask)
            elif binl[jbin]>binu[jbin]: # if the bin goes accross 180
                in_bin = (((R1>=binl[jbin]) | (R1<binu[jbin])) & QC & outlier_mask)
            
            if sum(in_bin)>min_in_bin:
                binStd = np.nanstd(A1[in_bin==1])
                binMedian = np.nanmedian(A1[in_bin==1])
                binMean = np.nanmean(A1[in_bin==1])
                binMean = np.average(A1[in_bin==1], weights=Weights_a[in_bin==1])
                binSum = sum(in_bin)
                binSum=len(np.unique(S0[in_bin]))
                
                

                if len(Weights_a)==np.sum(Weights_a):
                    binErr = binStd/np.sqrt(binSum)
                else:
                    binErr = 1/np.sqrt(np.sum(Weights_a[in_bin==1]))
                    # https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
                    
                    
                if BOOTSTRAP:
                    N_iter = 100
                    binMeans = np.arange(0,N_iter)*np.NaN
                    df = pd.DataFrame(data=np.transpose([S0[in_bin==1],A1[in_bin==1],np.square(Weights_a[in_bin==1])]))
                    df=(df.groupby(df[df.columns[0]]).mean())
                    for j in np.arange(0,N_iter):
                        #A1_rs, weights_rs = resample( A1[in_bin==1],Weights_a[in_bin==1])
                        
                        A1_rs, weights_rs = resample( df[df.columns[0]].values,np.sqrt(df[df.columns[1]].values), random_state=j+42)

                        
                        binMeans[j]=np.average(A1_rs, weights=weights_rs)
                        binErr = np.std(binMeans)
                
                    
            else:
                binMean=np.nan; binMedian=np.nan; binStd=np.nan; binSum=sum(in_bin); binErr=np.nan; binErr2=np.nan;
                binSum=len(np.unique(S0[in_bin]))
        else:
            binMean=np.nan; binMedian=np.nan; binStd=np.nan; binSum=sum(in_bin); binErr=np.nan; binErr2=np.nan;
            binSum=len(np.unique(S0[in_bin]))
            
        
        if binl[jbin]<binu[jbin]:
            bin_data['x_mean'].append(np.nanmean(R1[in_bin==1]))
        elif binl[jbin]>binu[jbin]: 
            # if the bin goes accross 180
            # need to bring lower part on top side, average and then ang180 again
            r1_Leq180 = R1[(((R1>=binl[jbin]) ) & QC & outlier_mask)]
            r1_Geq180 = R1[(((R1<binu[jbin]) ) & QC & outlier_mask)]
            bin_data['x_mean'].append(wvc.ang180  ( np.nanmean( np.concatenate([r1_Leq180, (r1_Geq180+360) ]) ) ) )
            
        bin_data['y_mean'].append(binMean)
        bin_data['y_median'].append(binMedian)
        bin_data['y_std'].append(binStd)
        bin_data['samples'].append(binSum)
        bin_data['y_err'].append(binErr)


    bin_data = pd.DataFrame(bin_data)
    bin_data['x_mean']=wvc.ang180(bin_data['x_mean'])
    bin_data = bin_data.sort_values('x_mean') # sort along R1 direction
    return bin_data, outlier_mask, outlier_count


def tryoshnikov_afc_unique(t,d,s,D,S,S0, QC, high_res=False,find_IQR_outlier=False, BOOTSTRAP=False, Weights_a=[], Weights_d=[]):
    """
        Function to run the flow distortion anlaysis for the ACE data acknowledging that the ERA-5 5-minute observations are not statistically independent
        with Options to 
        - identify and exclude outliers
        - calculate weighted averages (by specifying the weights of the observations)
        - apply bootstrapping to estimate the uncertainty of the mean
        - use disjunct or overlapping wind direction intervals

        :param t: date time index
        :param d: observed relative wind direction
        :param s: observed relative wind speed
        :param D: expected relative wind direction (calculated from the freestreem reference)
        :param S: expected relative wind speed (calculated from the freestreem reference)
        :param S0: variable used to identify if the observations are truely unique per wind direction sector (e.g. use ERA-5 U10N)
        :param QC: quality flag (only data [QC==True] are used)
        :high_res: True/False: if True, the wind direction bins are overlapping, if false, the wind direction bins are disjunct
        :param find_IQR_outlier: True/False, if True, data for each bin are screened for outlieres using outliers_iqr_noise(A1,NOISE)
        :param NOISE: noise level of A1 used in outliers_iqr_noise(A1,NOISE) to avoid rejection of digital noise
        :param BOOTSTRAP: True/Flase, if True the uncertainty is estimated via bootstrapping, if False, normal standard deviations are calcualted
        :param Weights_a: (OPTIONAL) Weights for s/S to be used in the averaging
        :param Weights_d: (OPTIONAL) Weights for d-D to be used in the averaging
        :returns: radqc_: Dataframe containing R=ang180(d),A=s/S,D=d-D,QC=QC,outliers={outliers identified via IQR filter} 
        :returns: afc_: Dataframe containing the bin mean, median and std values of R, A, D,
    """
    # run ratio and dirdiff analysis
    # optional filter based on IQR
    # winddirection bins are hardset inside this function, they need to be adapted to the dataset.
    min_in_bin = 7#12

    R = wvc.ang180(d)
    A = s/S
    dD = wvc.ang180(wvc.ang180(d)-wvc.ang180(D))
    if high_res:
        binl = np.concatenate( [np.arange(-180,-130,2), np.arange(-130,180,1)] ) # was at -140
        binu = np.concatenate( [np.arange(-180,-130,2)+10, np.arange(-130,180,1)+5] )
    else:
        binl = np.concatenate( [np.arange(-180,-130,10), np.arange(-130,185,5)] ) # was at -140
        binu = np.concatenate( [np.arange(-180,-130,10)+10, np.arange(-130,185,5)+5] )

    outliers = np.zeros_like(QC)

    bin_data_A, outlier_mask_A, outlier_count_A = Rdir_bin_intervals(R,A,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.05, Weights_a=Weights_a, BOOTSTRAP=False)
    bin_data_D, outlier_mask_D, outlier_count_D = Rdir_bin_intervals(R,dD,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.1, Weights_a=Weights_d, BOOTSTRAP=False)

    outliers = outliers | (outlier_count_A | outlier_count_A) # collect all outliers
    if find_IQR_outlier:
        QC = (QC & outlier_mask_A & outlier_mask_D)
    if 1: # always run a second time
        bin_data_A, outlier_mask_A, outlier_count_A = Rdir_bin_intervals_unique(R,A,S0,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.05, Weights_a=Weights_a, BOOTSTRAP=BOOTSTRAP)
        bin_data_D, outlier_mask_D, outlier_count_D = Rdir_bin_intervals_unique(R,dD,S0,QC,binl,binu,min_in_bin=min_in_bin,find_IQR_outlier=find_IQR_outlier,NOISE=0.1, Weights_a=Weights_d, BOOTSTRAP=BOOTSTRAP)
        QC = QC & (outlier_mask_A & outlier_mask_D)
        outliers = outliers | (outlier_count_A | outlier_count_A) # collect all outliers


    afc_ = pd.DataFrame({'R':bin_data_A['x_mean'],'A':bin_data_A['y_mean'],'D':bin_data_D['y_mean'], 
                         'A_median':bin_data_A['y_median'],'D_median':bin_data_D['y_median'],
                         'A_err':(bin_data_A['y_err']), 'D_err':(bin_data_D['y_err']), 'samples':bin_data_A['samples']})
    # afc data
    radqc_ = pd.DataFrame({'R':R,'A':A,'D':dD,'QC':QC,'outliers':outliers})
    

    return radqc_, afc_


def zeta_zu_zL_limited(LMO,zu):
    """
        Function to limit the ratio of z/L to 5
        This is taken from  ‘Part IV : Physical Processes’. In IFS Documentation CY45R1. IFS Documentation 4. ECMWF, 2018. https://www.ecmwf.int/node/18714.
        Section 3.2. THE SURFACE LAYER
        "... In extremely stable situations, i.e. for very small positive L, the ratio z/L is large, resulting in unrealistic
        profile shapes with standard stability functions. Therefore the ratio z/L is limited to 5 by defining a
        height h such that h=L = 5. If z < h, then the profile functions described above, are used up to z = h
        and the profiles are assumed to be uniform above that. This modification of the profiles for exceptionally
        stable situations (no wind) is applied to the surface transfer formulation as well as to the interpolation
        for post-processing."
        
        :param LMO:
        :param zu:
        
        :returns: zeta:
        :returns: zu:
    """
    zeta_limited = zu/LMO
    zu_limited = np.ones_like(zeta_limited)*zu
    zu_limited[zeta_limited>5]=5*LMO[zeta_limited>5]
    zeta_limited[zeta_limited>5]=5
    
    return zeta_limited, zu_limited

def expected_relative_wind(era5, wind_m, Zanemometer):
    """
        Function to calculate the expected relative wind speed (WSR) and expected relative wind direction (WDR) at the anemometer positions from the ERA-5 reference wind vector and the ships velocity and heading

        :param era5: data fram with cloumns ['ustar', 'WS10', 'u10', 'v10', 'LMO']
        :param wind_m: data fram with cloumns ['velEast', 'velNorth', 'HEADING']
        :param Zanemometer: height of the anemometer above sea level [meter]
         
        :returns: era5: Dataframe era5 with additional columns ['Urel', 'Vrel', 'WSR', 'WDR', 'WDIR', 'WS30']
    """
    # calculate U30 (z=31.5)
    #zeta30 = Zanemometer/era5['LMO']
    #z30 = np.ones_like(era5['LMO'])*Zanemometer
    #if 1:
    #    z30[zeta30>5]=5*era5['LMO'][zeta30>5]
    #    zeta30[zeta30>5]=5
        
    zeta30, z30 = zeta_zu_zL_limited(era5['LMO'],Zanemometer)
    era5['WS30'] = aceairsea.coare_u2ustar (era5['ustar'], input_string='ustar2u', coare_version='coare3.5', TairC=(era5.t2m-273.15), z=z30, zeta=zeta30)

    
    era5 = era5.assign( Urel = (era5.u10*era5.WS30/era5.WS10 - wind_m.velEast));  # relative wind speed in earth frame
    era5 = era5.assign( Vrel = (era5.v10*era5.WS30/era5.WS10 - wind_m.velNorth)); # relative wind speed in earth frame
    era5 = era5.assign( WSR = (np.sqrt(era5.Urel*era5.Urel+era5.Vrel*era5.Vrel)));
    era5 = era5.assign( WDR = ((-wind_m.HEADING + 270 - np.rad2deg(np.arctan2(era5.Vrel, era5.Urel)) )% 360 ));    
    era5 = era5.assign( WDIR=((270-np.rad2deg(np.arctan2(era5.v10,era5.u10))) % 360) ) # coming from direction

    return era5

def flag_4DVAR_affected_data(UBXH3_assimilated, wind_m):
    """
        Function to read the list of assimilation events and derive a flag to exclude potentially affected data

        :returns: a dataframe with a flag that can be used to identify observations that may be affected by 4DVAR assimilation of wind data that was reported from the Akademik Tryoshnikov under the sation id UBXH3
    """
    # returns data frame with unified time stamp and column '4DVAR' = 1 if data is affected 
    # assume 9:00 to 21:00 and 21:00 to 9:00 windows are affected if one reading is within
    UBXH3_assimilated['4DVAR']=1
    UBXH3_4DVAR = pd.DataFrame(index=wind_m.index, columns=[])
    UBXH3_4DVAR = UBXH3_4DVAR.merge( UBXH3_assimilated.append(pd.DataFrame(index=(UBXH3_assimilated.index+pd.to_timedelta(12, unit='h')), columns=UBXH3_assimilated.columns ).iloc[-1]).resample('12H', loffset = datetime.timedelta(hours=(-3))).mean().resample('5T', loffset = datetime.timedelta(seconds=(5*60/2) ) ).mean().interpolate(limit=12*12, limit_direction='forward')[['4DVAR']], left_index=True,right_index=True, how='left' )
    if 0:
        fig = plt.figure()
        fig.set_size_inches(20,8)
        plt.plot(UBXH3_assimilated['4DVAR'],'o')
        plt.plot(UBXH3_assimilated.resample('12H', loffset = datetime.timedelta(hours=(-3))).mean().resample('1H').mean().interpolate(limit=12, limit_direction='forward')['4DVAR'], '.' )
        plt.plot(UBXH3_4DVAR['4DVAR'], 'k')
        #plt.xlim([pd.to_datetime('2017-01-26 08:00:00'), pd.to_datetime('2017-01-28 18:00:00')])
        #plt.xlim([pd.to_datetime('2017-02-07 00:00:00'), pd.to_datetime('2017-02-08 18:00:00')])
        #plt.xlim([pd.to_datetime('2017-02-18 00:00:00'), pd.to_datetime('2017-02-19 00:00:00')])
    
    return UBXH3_4DVAR




if __name__ == "__main__":
    print('running the air flow distortion estimation and correction for the ACE data set')
    from pathlib import Path

    # import local functionalities
    from pyantarctica import aceairsea as aceairsea 
    from pyantarctica import windvectorcoordinates as wvc
    import read_ace_data as read_ace_data
    
    
    FOLD_out = './data/flow_distortion_bias/'
    AFC_BASE_FILE_NAME = 'flow_distortion_bias_sensor'
    afc_correction_factor_files = str(Path(FOLD_out, AFC_BASE_FILE_NAME))
    plot_folder = './plots/'
    BOOTSTRAP = True # flag to calculate bin error form weighter mean formular of via bootstrap
    HIGH_RES = True # change this to true for better resolution of the correction factors (rejects 2% more outliers)
    
    Zanemometer = 31.5 # estimated height of the anemometer above mean sea level [meter]
    
    # create the output folders if not existent
    Path(FOLD_out).mkdir(parents=True, exist_ok=True)
    Path('./data/wind_data_corrected_fivemin/').mkdir(parents=True, exist_ok=True)
    Path('./data/wind_data_corrected_onemin/').mkdir(parents=True, exist_ok=True)
    Path('./data/wind_data_uncorrected_fivemin/').mkdir(parents=True, exist_ok=True)
    Path('./data/wind_data_uncorrected_onemin/').mkdir(parents=True, exist_ok=True)
    Path('./data/wind_data_corrected_combined_fivemin/').mkdir(parents=True, exist_ok=True)
        

    Path(plot_folder).mkdir(parents=True, exist_ok=True)

    ###################################
    # read all required data
    ###################################
    
    # wind/gps-velocity data at 1 minute
    wind_m = read_ace_data.wind_merge_gps_afc_option(afc_correction_factor_files=[])
    
    if 1:
        # resample to 1 minutes and save the result
        wind_1min = wvc.resample_wind_data(wind_m, Nmin=1, interval_center='odd', lon_flip_tollerance=0.0005)
        wind_1min.latitude = wind_1min.latitude.interpolate()
        wind_1min.longitude = wind_1min.longitude.interpolate()

        wind_m_CF_stbd = wind_1min[['latitude','longitude', 'WDR1', 'WSR1','WD1','WS1','uR1', 'vR1', 'u1', 'v1']].copy()
        wind_m_CF_port = wind_1min[['latitude','longitude', 'WDR2', 'WSR2','WD2','WS2','uR2', 'vR2', 'u2', 'v2']].copy()

        wind_m_CF_stbd = wind_m_CF_stbd.rename(columns={'WDR1':'wind_from_direction_relative_to_platform', 
                                                        'WSR1':'wind_speed_relative_to_platform',
                                                        'WD1' : 'wind_from_direction',
                                                        'WS1' : 'wind_speed',
                                                        'uR1' : 'bowward_relative_wind',
                                                        'vR1' : 'portward_relative_wind',
                                                        'u1' : 'eastward_wind',
                                                        'v1' : 'northward_wind',
                                                       })
        wind_m_CF_port = wind_m_CF_port.rename(columns={'WDR2':'wind_from_direction_relative_to_platform', 
                                                        'WSR2':'wind_speed_relative_to_platform',
                                                        'WD2' : 'wind_from_direction',
                                                        'WS2' : 'wind_speed',
                                                        'uR2' : 'bowward_relative_wind',
                                                        'vR2' : 'portward_relative_wind',
                                                        'u2' : 'eastward_wind',
                                                        'v2' : 'northward_wind',
                                                       })

        wind_m_CF_stbd.to_csv('./data/wind_data_uncorrected_onemin/wind-observations-stbd-uncorrected-1min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
        wind_m_CF_port.to_csv('./data/wind_data_uncorrected_onemin/wind-observations-port-uncorrected-1min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')


    
    
    # resample to 5 minutes
    wind_m = wvc.resample_wind_data(wind_m, Nmin=5, interval_center='odd', lon_flip_tollerance=0.01)
    # interpolate the 5min lat lon this does a good job.
    wind_m.latitude = wind_m.latitude.interpolate()
    wind_m.longitude = wind_m.longitude.interpolate()
    
    if 1:
        # save the 5-minute data
        wind_m_CF_stbd = wind_m[['latitude','longitude', 'WDR1', 'WSR1','WD1','WS1','uR1', 'vR1', 'u1', 'v1']].copy()
        wind_m_CF_port = wind_m[['latitude','longitude', 'WDR2', 'WSR2','WD2','WS2','uR2', 'vR2', 'u2', 'v2']].copy()

        wind_m_CF_stbd = wind_m_CF_stbd.rename(columns={'WDR1':'wind_from_direction_relative_to_platform', 
                                                        'WSR1':'wind_speed_relative_to_platform',
                                                        'WD1' : 'wind_from_direction',
                                                        'WS1' : 'wind_speed',
                                                        'uR1' : 'bowward_relative_wind',
                                                        'vR1' : 'portward_relative_wind',
                                                        'u1' : 'eastward_wind',
                                                        'v1' : 'northward_wind',
                                                       })
        wind_m_CF_port = wind_m_CF_port.rename(columns={'WDR2':'wind_from_direction_relative_to_platform', 
                                                        'WSR2':'wind_speed_relative_to_platform',
                                                        'WD2' : 'wind_from_direction',
                                                        'WS2' : 'wind_speed',
                                                        'uR2' : 'bowward_relative_wind',
                                                        'vR2' : 'portward_relative_wind',
                                                        'u2' : 'eastward_wind',
                                                        'v2' : 'northward_wind',
                                                       })
        wind_m_CF_stbd.to_csv('./data/wind_data_uncorrected_fivemin/wind-observations-stbd-uncorrected-5min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
        wind_m_CF_port.to_csv('./data/wind_data_uncorrected_fivemin/wind-observations-port-uncorrected-5min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
    
    era5 = read_ace_data.read_era5_data()
    dist2land = read_ace_data.read_distance2land()
    UBXH3_assimilated = read_ace_data.read_assimilation_list()
    
    
    UBXH3_4DVAR = flag_4DVAR_affected_data(UBXH3_assimilated, wind_m)

    if 1:
        era5 = pd.merge(era5, wind_m[['HEADING']], left_index=True,right_index=True, how='right')
        era5.drop(columns=['HEADING'], inplace=True)
    
        dist2land = pd.merge(dist2land, wind_m[['HEADING']], left_index=True,right_index=True, how='right')
        dist2land.drop(columns=['HEADING'], inplace=True)

    # calculate the expected relative wind speed and direction at anemometer position and add these fileds to the era5 Dataframe
    era5 = expected_relative_wind(era5, wind_m, Zanemometer)
        
    ###################################
    # estimate expected error in WSR_Model
    ###################################
    
    WSR_err, WDR_err = wvc.WSRWDR_uncertainy(era5.WS30,era5.WDIR,wind_m.HEADING,wind_m.velEast,wind_m.velNorth,a_WSPD=1.2,d_WDIR=10)
    a1_err = np.sqrt(np.square(0.01)+np.square(WSR_err*wind_m.WSR1/era5.WSR/era5.WSR) ) # the (absolute) error in a1 estimate is given by
    Weights_a1=(1/np.square(a1_err))
    a2_err = np.sqrt(np.square(0.01)+np.square(WSR_err*wind_m.WSR2/era5.WSR/era5.WSR) ) # the (absolute) error in a1 estimate is given by
    Weights_a2=(1/np.square(a2_err))

    d_err = np.sqrt(np.square(10)+np.square(WDR_err) ) # the (absolute) error in d estimate is given by
    Weights_d1=(1/np.square(d_err))
    Weights_d2=(1/np.square(d_err))
    
    # filter data based on stability of the ships speed and heading within the 5-minute intervals
    WDIR_MAX_DIFF = 15 # maximum variation between the 1min wind directions
    QC_WDIR_DIFF = (np.max([wind_m.WDR1_DIFF, wind_m.WDR2_DIFF], axis=0) < WDIR_MAX_DIFF)
    QC_WDIR1_DIFF = ((wind_m.WDR1_DIFF < WDIR_MAX_DIFF) & (wind_m.WSR1>2) )
    QC_WDIR2_DIFF = ((wind_m.WDR2_DIFF < WDIR_MAX_DIFF) & (wind_m.WSR2>2) )
    QC0_SOG_HDG = ( ( (wind_m.HEADING_DIFF<10) & (wind_m.SOG>2) ) | ( (wind_m.HEADING_DIFF<10) & (wind_m.SOG<2) ) )  & ( (wind_m.SOG_DIFF<1) )

    print('fraction of samples remaining after Wind Direction Variability check row 1 = stbd, row 2 = port')
    print('[fraction, remaining, total]')
    print([np.sum((QC_WDIR_DIFF==True) & (wind_m.WSR1>-1) )/np.sum(wind_m.WSR1>-1), np.sum((QC_WDIR_DIFF==True)), np.sum(wind_m.WSR1>-1)])
    print([np.sum((QC_WDIR_DIFF==True) & (wind_m.WSR2>-1) )/np.sum(wind_m.WSR2>-1), np.sum((QC_WDIR_DIFF==True)), np.sum(wind_m.WSR2>-1)])
    
    print("Calculating the relative wind speed and direction ratio of sensor 1 and 2 (starboard and port sensor) ...")
    radqc_, afc_ = tryoshnikov_afc(wind_m.index,wind_m.WDR1,wind_m.WSR1,wind_m.WDR2,wind_m.WSR2, QC=(QC_WDIR_DIFF ), high_res=False,find_IQR_outlier=True, BOOTSTRAP=BOOTSTRAP)
    print("... done!")
        
    QC_ECMWF = (era5.WSR>2) & (era5.LSM==0) & (era5.SIF==0) & (Zanemometer/era5.LMO>-1.5) &  (Zanemometer/era5.LMO<.25) #10%correction
    print('QC_ECMWF==1: '+str(np.sum(QC_ECMWF)) )
    QC_ECMWF = QC_ECMWF & (dist2land.distance>50000) # could also use 25000
    print('QC_ECMWF==1: '+str(np.sum(QC_ECMWF)) )
    QC1 = (QC_WDIR1_DIFF & QC0_SOG_HDG & QC_ECMWF & (radqc_.outliers==0) & np.isnan(UBXH3_4DVAR['4DVAR']) ) # 
    QC2 = (QC_WDIR2_DIFF & QC0_SOG_HDG & QC_ECMWF & (radqc_.outliers==0) & np.isnan(UBXH3_4DVAR['4DVAR']) )
    print('port: QC passed by '+str(np.sum(QC1))+' of '+ str(np.sum(wind_m.WS1>-1)) + ' samples; fraction of data used '+  str(np.sum(QC1)/np.sum(wind_m.WS1>-1)) )
    print('stbd: QC passed by '+str(np.sum(QC2))+' of '+ str(np.sum(wind_m.WS1>-1)) + ' samples; fraction of data used '+  str(np.sum(QC2)/np.sum(wind_m.WS2>-1)) )

    
    ###################################
    # calculate the correction factors
    ###################################
    
    print("Calculating the flowdistortion bias of sensor 1 (starboard)")
    radqc_s1, afc_s1 = tryoshnikov_afc_unique(wind_m.index,wind_m.WDR1,wind_m.WSR1,era5.WDR,era5.WSR,S0=era5.WS10N, QC=QC1, high_res=HIGH_RES,find_IQR_outlier=True, BOOTSTRAP=BOOTSTRAP, Weights_a=Weights_a1, Weights_d=Weights_d1  )
    print("... done!")
    
    print("Calculating the flowdistortion bias of sensor 2 (port)")
    radqc_s2, afc_s2 = tryoshnikov_afc_unique(wind_m.index,wind_m.WDR2,wind_m.WSR2,era5.WDR,era5.WSR,S0=era5.WS10N, QC=QC2, high_res=HIGH_RES,find_IQR_outlier=True, BOOTSTRAP=BOOTSTRAP, Weights_a=Weights_a2, Weights_d=Weights_d2  )
    print("... done!")
    
    afc_s1 = afc_s1.rename(columns={'R':'wind_from_direction_relative_to_platform', 
                            'A': 'mean_of_wind_speed_bias',
                            'D': 'mean_of_wind_direction_bias',
                            'A_median': 'median_of_wind_speed_bias',
                            'D_median': 'median_of_wind_direction_bias',
                            'A_err': 'uncertainty_of_wind_speed_bias',
                            'D_err': 'uncertainty_of_wind_direction_bias',
                            'samples': 'number_of_samples'})

    afc_s2 = afc_s2.rename(columns={'R':'wind_from_direction_relative_to_platform', 
                        'A': 'mean_of_wind_speed_bias',
                        'D': 'mean_of_wind_direction_bias',
                        'A_median': 'median_of_wind_speed_bias',
                        'D_median': 'median_of_wind_direction_bias',
                        'A_err': 'uncertainty_of_wind_speed_bias',
                        'D_err': 'uncertainty_of_wind_direction_bias',
                        'samples': 'number_of_samples'})
    
    afc_s1.to_csv( Path( FOLD_out + AFC_BASE_FILE_NAME + '1.csv') , index=False , float_format='%.4f')
    afc_s2.to_csv( Path( FOLD_out + AFC_BASE_FILE_NAME + '2.csv') , index=False , float_format='%.4f')
    
    ###################################
    # calculate the corrected wind speed and direction
    ###################################
    
    wind_c = read_ace_data.wind_merge_gps_afc_option(afc_correction_factor_files=afc_correction_factor_files)
    
    if 1:
        # resample to 1 minutes overwrite the same dataframe to safe space
        wind_1min = wvc.resample_wind_data(wind_c, Nmin=1, interval_center='odd', lon_flip_tollerance=0.0005)
        wind_1min.latitude = wind_1min.latitude.interpolate()
        wind_1min.longitude = wind_1min.longitude.interpolate()

        wind_c_CF_stbd = wind_1min[['latitude','longitude', 'WDR1', 'WSR1','WD1','WS1','uR1', 'vR1', 'u1', 'v1']].copy()
        wind_c_CF_port = wind_1min[['latitude','longitude', 'WDR2', 'WSR2','WD2','WS2','uR2', 'vR2', 'u2', 'v2']].copy()

        wind_c_CF_stbd = wind_c_CF_stbd.rename(columns={'WDR1':'wind_from_direction_relative_to_platform', 
                                                        'WSR1':'wind_speed_relative_to_platform',
                                                        'WD1' : 'wind_from_direction',
                                                        'WS1' : 'wind_speed',
                                                        'uR1' : 'bowward_relative_wind',
                                                        'vR1' : 'portward_relative_wind',
                                                        'u1' : 'eastward_wind',
                                                        'v1' : 'northward_wind',
                                                       })
        wind_c_CF_port = wind_c_CF_port.rename(columns={'WDR2':'wind_from_direction_relative_to_platform', 
                                                        'WSR2':'wind_speed_relative_to_platform',
                                                        'WD2' : 'wind_from_direction',
                                                        'WS2' : 'wind_speed',
                                                        'uR2' : 'bowward_relative_wind',
                                                        'vR2' : 'portward_relative_wind',
                                                        'u2' : 'eastward_wind',
                                                        'v2' : 'northward_wind',
                                                       })

        wind_c_CF_stbd.to_csv('./data/wind_data_corrected_onemin/wind-observations-stbd-corrected-1min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
        wind_c_CF_port.to_csv('./data/wind_data_corrected_onemin/wind-observations-port-corrected-1min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')

    
    
    wind_c = wvc.resample_wind_data(wind_c, Nmin=5,interval_center='odd', lon_flip_tollerance=0.01)

    if 1:
        wind_c_CF_stbd = wind_c[['latitude','longitude', 'WDR1', 'WSR1','WD1','WS1','uR1', 'vR1', 'u1', 'v1']].copy()
        wind_c_CF_port = wind_c[['latitude','longitude', 'WDR2', 'WSR2','WD2','WS2','uR2', 'vR2', 'u2', 'v2']].copy()

        wind_c_CF_stbd = wind_c_CF_stbd.rename(columns={'WDR1':'wind_from_direction_relative_to_platform', 
                                                        'WSR1':'wind_speed_relative_to_platform',
                                                        'WD1' : 'wind_from_direction',
                                                        'WS1' : 'wind_speed',
                                                        'uR1' : 'bowward_relative_wind',
                                                        'vR1' : 'portward_relative_wind',
                                                        'u1' : 'eastward_wind',
                                                        'v1' : 'northward_wind',
                                                       })
        wind_c_CF_port = wind_c_CF_port.rename(columns={'WDR2':'wind_from_direction_relative_to_platform', 
                                                        'WSR2':'wind_speed_relative_to_platform',
                                                        'WD2' : 'wind_from_direction',
                                                        'WS2' : 'wind_speed',
                                                        'uR2' : 'bowward_relative_wind',
                                                        'vR2' : 'portward_relative_wind',
                                                        'u2' : 'eastward_wind',
                                                        'v2' : 'northward_wind',
                                                       })
        wind_c_CF_stbd.to_csv('./data/wind_data_corrected_fivemin/wind-observations-stbd-corrected-5min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
        wind_c_CF_port.to_csv('./data/wind_data_corrected_fivemin/wind-observations-port-corrected-5min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')
    
    
    # take the average of sensor 1 and sensor 2
    u = np.nanmean([wind_c.u1, wind_c.u2],axis=0);
    v = np.nanmean([wind_c.v1, wind_c.v2],axis=0);
    uR = np.nanmean([wind_c.uR1, wind_c.uR2],axis=0);
    vR = np.nanmean([wind_c.vR1, wind_c.vR2],axis=0);
    uz = np.sqrt(np.square(u)+np.square(v))
    dir10 = (270 - np.rad2deg( np.arctan2( v, u ) ) ) % 360

    # estimate u10N using the stability frm ERA-5
    TairC=era5.T2M
    
    zeta30, z30 = zeta_zu_zL_limited(era5['LMO']*np.square(uz/era5['WS30']),Zanemometer) # scale LMO with u10ratio^2  # this adjustment avoids overshooting corrections
    ustar = aceairsea.coare_u2ustar(uz, input_string='u2ustar', coare_version='coare3.5', TairC=TairC, z=z30, zeta=zeta30)
    u10N = aceairsea.coare_u2ustar(ustar, input_string='ustar2u', coare_version='coare3.5', TairC=TairC, z=10, zeta=0)

    wind_c_CF = wind_c_CF_stbd.copy()
    wind_c_CF['wind_from_direction_relative_to_platform'] = (180-np.rad2deg(np.arctan2(vR,uR) ) )%360
    wind_c_CF['wind_speed_relative_to_platform'] = np.sqrt(np.square(uR)+np.square(vR))
    wind_c_CF['wind_from_direction'] = dir10
    wind_c_CF['wind_speed'] = uz
    wind_c_CF['bowward_relative_wind'] = uR
    wind_c_CF['portward_relative_wind'] = vR
    wind_c_CF['eastward_wind'] = u
    wind_c_CF['northward_wind'] = v
    wind_c_CF['speed_of_10m_neutral_wind'] = u10N
    wind_c_CF['eastward_component_of_10m_neutral_wind'] = u10N/uz*u
    wind_c_CF['northward_component_of_10m_neutral_wind'] = u10N/uz*v

    wind_c_CF.to_csv('./data/wind_data_corrected_combined_fivemin/wind-observations-port-stbd-corrected-combined-5min-legs0-4.csv',date_format="%Y-%m-%dT%H:%M:%S+00:00",na_rep="NaN", float_format='%.4f')

    
    




