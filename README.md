# wind-speed-correction

Analysis of the flowdistortion bias in the ship borne wind speed measurements from the [Antarctic Circumnavigation Experiment (ACE)](https://doi.org/10.5281/zenodo.1443511). We use ERA-5 reanalysis data to quantify the flowdistortion bias in the measured relative wind speed and direction as a function of the relative wind direction. 

Please refer to and cite Landwehr et al. (2019; DOI: https://doi.org/10.5194/amt-2019-366) when using this code or modified versions of it.

## Using the code

An example of how these scripts can be used can be found in a [Gitlab repository](https://renkulab.io/gitlab/ACE-ASAID/wind-speed-correction) which also contains datasets that were used and created by this code.

### Requirements

In order to run the code, the packages in ```requirements.txt``` should be installed:

```pip3 install -r requirements.txt```

Note that [pyantarctica](https://github.com/Swiss-Polar-Institute/pyantarctica) requires an older version of pandas, but this code will run with the newer version of pandas=1.0.3.

### Running the code

Amend paths to input and output dataset files.

Then execute with python3:

```python3 ace_flowdistortion.py```

### Input datasets

This repository can be used with the following input datasets.

Summary raw wind data from the Southern Ocean collected on board the Antarctic Circumnavigation Expedition (ACE) during the austral summer of 2016/2017. https://doi.org/10.5281/zenodo.3801718

One-minute average cruise track and ship velocity of the Antarctic Circumnavigation Expedition (ACE) undertaken during the austral summer of 2016/2017. https://doi.org/10.5281/zenodo.3772895

ERA-5 reanalysis results interpolated onto the five-minute average cruise track of the Antarctic Circumnavigation Expedition (ACE) during the austral summer of 2016/2017. https://doi.org/10.5281/zenodo.3831980

Distance to the nearest land/coastline (including small subantarctic islands) for the five-minute average cruise track of the Antarctic Circumnavigation Expedition (ACE) during the austral summer of 2016/2017. https://doi.org/10.5281/zenodo.3832045

77040_UBXH3_2017012618_2017021818.txt provides a list of all the times that near surface wind observations were assimilated into ERA5 from station UBXH3 (the R/V Akademik Tryoshnikov used during ACE), during the period November 2016 to May 2017. The columns in the file are stationId, date (yyyymmdd; UTC), time (hh; UTC), latitude (degree_north), longitude (degree_east).

All of the above datasets are available with a Creative Commons Attribution 4.0 International License (CC BY 4.0) whose full text can be found at https://creativecommons.org/licenses/by/4.0/

### Output datasets

Examples of the output datasets can be found in the Gitlab repositories or in the following publications: 
TODO
