import numpy as np
import pandas as pd
from datetime import datetime, timedelta

def add_legs_index(df, **kwargs):

    """
        Add a leg identifier (defaults to [1,2,3]) to a datetime indexed dataframe of ACE measurements

        :param df: datetime indexed dataframe to which add the leg index
        :param kwargs: contains options to override defaults:

        |   leg_dates : array indicating beginning and end (columns) of time period belonging to a leg / subleg (number of rows equals number of legs / sublegs)
        |   codes : list indicating what code to use to denote legs. Defaults to integers in range [1,2,3,...] with 0 indicating out-of-leg datapoints

        :returns: dataframe with a new column "leg" indicating to which leg the datapoint belongs to

    """

    if "leg_dates" not in kwargs:
        # leg_dates = [['2016-12-20', '2017-01-21'], # leg 1
        #             ['2017-01-22', '2017-02-25'],  # leg 2
        #             ['2017-02-26', '2017-03-19']]  # leg 3
        # updates leg dates with hourly resolution
        leg_dates = [
            ["2016-11-17", "2016-12-19"],  # leg 0
            ["2016-12-20 16:00", "2017-01-18 22:00"],  # leg 1
            ["2017-01-22 11:00", "2017-02-22 12:00"],  # leg 2
            [
                "2017-02-26 01:00",
                "2017-03-19 09:00",
            ],  # leg 3, ship at full speed @'2017-02-26 02:00', in the vicinity at '2017-03-18 14:00'
            ["2017-03-22 19:00", "2017-04-11 16:00"],  # leg 4
        ]

    else:
        leg_dates = kwargs["leg_dates"]

    # merz_glacier = ['2017-01-29', '2017-01-31']

    if "codes" not in kwargs:
        codes = np.arange(1, 1 + len(leg_dates))
        codes = [0, 1, 2, 3, 4]  # SL
    else:
        codes = kwargs["codes"]

    """Add a column to the datatable specifying the cruise leg"""
    assert len(codes) == len(
        leg_dates
    ), "To each date interval must correspond only one code"

    if "leg" in df:
        print("leg column already there")
        return df

    dd = pd.Series(data=np.zeros((len(df.index))) * np.nan, index=df.index, name="leg")

    c = 0
    while c < len(codes):
        dd.loc[leg_dates[c][0] : leg_dates[c][1]] = codes[c]
        c += 1

    if "merz_glacier" in kwargs:
        dd.loc[merz_glacier[0] : merz_glacier[1]] = -1

    #  dd.loc[dd['leg'].isnull()] = 0
    df = df.assign(leg=dd)

    return df


##############################################################################################################
def match2series(ts, ts2match):
    """
        Function to crop/append to the series ts in order to have the same number of samples as ts2match
        REQUIRES the two indicees to be same for ts2match and ts !!!
        TODO add a warning if this is not the case!!!
        this is just a wrapper around pandas merge (import pandas as pd)

        :param ts: datetime indexed dataframe to be modified
        :param ts2match: datetime indexed dataframe, of which the index shall be taken

        :Example:

            wind = match2series(wind,aerosols) # output is the wind matched to aerosols
    """
    ts2match = ts2match[ts2match.columns[0]].to_frame()
    ts2match.rename(index=str, columns={ts2match.columns[0]: "var2match"}, inplace=True)
    ts = pd.merge(ts, ts2match, left_index=True, right_index=True, how="right")
    ts = ts.drop(columns=["var2match"])
    return ts


##############################################################################################################


def resample_timeseries(
    ts,
    time_bin,
    how="mean",
    new_label_pos="c",
    new_label_parity="even",
    old_label_pos="c",
    old_resolution=0,
    COMMENTS=False,
):
    """
        Function to resample timeseries to multiple of minutes placing the inter val label where you like it
        Info on time series location in the inital ts can be used to ensure accurate binning
        Output time stamp label position left, right, centre and parity ('even' -> index=00:05:00, 'odd' -> index=00:02:30) can be choosen

        :param ts: datetime indexed dataframe to resample
        :param time_bin: integer aggregation time, in minutes
        :param how: string specifyin how to aggregate. Has to be compatible with df.resample('5T').aggregate(how)
        :param old_label_pos: string ('l'=left, 'r'=right, 'c'=center) position of the initial timestamp
        :param old_resolution: integer input time resolution in minutes used to correct input time stamp if (old_label_pos=='c')==False, set to 0 if unknown
        :param new_label_pos: string ('l'=left, 'r'=right, 'c'=center), define if timest_ label denotes left, right, center of new intervals
        :param new_label_parity: string ('even' -> index=00:05:00, 'odd' -> index=00:02:30), if time stamp will look like as when resample would be run ('even')
        :param COMMENTS: boolean (True->print some info about what is done)
        :returns: resampled dataframe with uninterrupted datetime index and NaNs where needed
        :Example:

            ts_mean = dataset.resample_timeseries(ts, 15, how='mean')
            # if you know that initial ts lable was on right of interval and resolution was 5min (e.g. aerosols)
            ts_mean = dataset.resample_timeseries(ts, 15, how='mean', old_label_pos='r', old_resolution=5)

    """

    # assume input time series has label on interval center, then
    # we can resample to new time series with label position in center of interval,
    # choose if you like the lable to look like 'even' -> index=00:05:00, 'odd' -> index=00:02:30
    if (new_label_pos == "c") & (new_label_parity == "even"):
        # put label on center, index=00:05:00
        ts_offset = timedelta(minutes=(time_bin / 2))
        rs_loffset = timedelta(minutes=0)
    elif (new_label_pos == "c") & (new_label_parity == "odd"):
        # put label on center, index=00:02:30
        ts_offset = timedelta(minutes=0)
        rs_loffset = timedelta(minutes=(time_bin / 2))
    elif (new_label_pos == "l") & (new_label_parity == "even"):
        # put label on left, index=00:05:00 (classic resample behaviour)
        ts_offset = timedelta(minutes=0)
        rs_loffset = timedelta(minutes=0)
    elif (new_label_pos == "l") & (new_label_parity == "odd"):
        # put label on left, index=00:02:30
        ts_offset = timedelta(minutes=-(time_bin / 2))
        rs_loffset = timedelta(minutes=+(time_bin / 2))
    elif (new_label_pos == "r") & (new_label_parity == "even"):
        # put label on right end of new resample interval, index=00:05:00
        ts_offset = timedelta(minutes=+time_bin)
        rs_loffset = timedelta(minutes=0)
    elif (new_label_pos == "r") & (new_label_parity == "odd"):
        # put label on right end of new resample interval, index=00:02:30
        ts_offset = timedelta(minutes=0)
        rs_loffset = timedelta(minutes=+time_bin)
    else:
        print('new_label_pos must be either "l","r", or "c"!')
        print('new_label_parity must be either "odd" or "even"!')
        return

    # now check if the old lable pos is not 'c' and add an offset to ts_offset to correct for this

    if (old_label_pos == "c") == False:
        if old_resolution > 0:
            # known old_resolution we can use it to calcualte the offset to add
            if old_label_pos == "l":
                ts_offset = ts_offset + timedelta(minutes=+(old_resolution / 2))
            elif old_label_pos == "r":
                ts_offset = ts_offset + timedelta(minutes=-(old_resolution / 2))
        else:
            tres_ = np.median(np.diff(ts.index.tolist())).total_seconds()
            tres_ = int(tres_)  # round to full second
            if COMMENTS == True:
                print(
                    "Inferring old time resolution to be "
                    + str(tres_)
                    + " seconds ("
                    + str(tres_ / 60)
                    + " minutes)"
                )
            if old_label_pos == "l":
                ts_offset = ts_offset + timedelta(seconds=+(tres_ / 2))
            elif old_label_pos == "r":
                ts_offset = ts_offset + timedelta(seconds=-(tres_ / 2))
            else:
                print('old_label_pos must be either "l","r", or "c"')
                return

    # fix the initial index if it is needed,
    ts_resampled = ts.copy()
    ts_resampled.index = ts_resampled.index + ts_offset
    # resample with desired offset (the loffset changes the lable after resample has acted on the time series)
    ts_resampled = ts_resampled.resample(
        str(time_bin) + "T", loffset=rs_loffset
    ).aggregate(how)
    return ts_resampled