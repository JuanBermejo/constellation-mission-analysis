# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:18:57 2021

@author: J.Álvarez & J.Bermejo
"""

#Basic libraries
import numpy as np
import pandas as pd
import os

# #For TLE fetching
# from satellite_tle import fetch_tle_from_celestrak

#For satellite propagation
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS

#Units, and coordinates conversion
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation


rel_path = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))


def get_az_sequence(TLE, GS, date, period):
    """
    This function computes the date of the accesses for the period
    demanded and the sequence of azimuth to return the rotor safely
    to 0.

    Parameters
    ----------
    date :  str
        Date in ISO format: "2021-01-01 12:00:00.000"
    period :  float
        Extension of time considered to detect accesses (days)

    Returns
    -------
    access_df :  pandas DataFrame
        Accesses data and azimuth sequence
    """

    #configuration parameters and variable inicialization
    dT_coarse = 60.
    dT_fine = 1.
    date0 = Time(date, format='iso')
    access_duration_estimated = 25/(24 * 60) #fraction of days
    access_list = []

    TLE_1 = TLE[0]
    TLE_2 = TLE[1]
    satellite = Satrec.twoline2rv(TLE_1, TLE_2)

    #Ground Station data

    lon_GS = GS[0]
    lat_GS = GS[1]
    height_GS = GS[2]
    min_elevation_GS = GS[3]

    #Coarse propagation
    r_coarse_teme, v_coarse_teme, time_coarse = propagate(satellite, date0, period, dT_coarse)

    #Compute azimuth and elevation
    _, el_coarse = teme2AzEl(r_coarse_teme, v_coarse_teme, time_coarse, lat_GS, lon_GS, height_GS)

    #Detect accesses
    idx_start, _ = detect_access(el_coarse, min_elevation_GS)

    ids = np.arange(len(idx_start))

    for idx in idx_start:
        row_dict = {}

        #Start of fine propagation 5 minutes before access start detected coarsely
        date_access = time_coarse[idx] - 5/(24 * 60)* u.day

        #Fine propagation for each access
        r_fine_teme, v_fine_teme, time_fine = propagate(satellite, date_access, access_duration_estimated, dT_fine)
        az_fine, el_fine = teme2AzEl(r_fine_teme, v_fine_teme, time_fine, lat_GS, lon_GS, height_GS)

        #index of start and end of access
        idx_start_fine, idx_end_fine = detect_access(el_fine, min_elevation_GS)

        start_date = (time_fine[idx_start_fine][0])
        end_date = (time_fine[idx_end_fine][0])
        duration = (end_date - start_date).value

        row_dict["Satellite"] = TLE_1[2:8]
        row_dict["Start Time [s]"] = (start_date - date0).value * 86400
        row_dict["Stop Time [s]"] = (end_date - date0).value * 86400
        row_dict["Duration [s]"] = duration * 86400
        # str(int(duration * 24 * 60)) + ":" + str(int(duration * 24 * 60  % 1 * 60))


        if any(el_fine > min_elevation_GS):
            access_list.append(row_dict)
        else:
            ids = ids[:-1]

    access_df = pd.DataFrame(data=access_list, index=ids)

    return access_df


def detect_access(elevation_array, min_elevation):
    """
    Detect start and end of access by checking elevation change and
    return the indexes where they happen.

    Parameters
    ----------
    elevation_array :  1D
        Date in ISO format: "2021-01-01 12:00:00.000"
    period :  float
        Extension of time considered to detect accesses (days)

    Returns
    -------
    idx_start :  1D-array
        Indexes of accesses start
    idx_start :  1D-array
        Indexes of accesses end
    """

    #add zero at the start and bottom of the array
    idx_ext = np.concatenate([[-1], elevation_array - min_elevation, [-1]])

    #detect start and end of an acces
    access_io = np.where(np.sign(idx_ext[:-1]) != np.sign(idx_ext[1:]))[0]

    #extract index of the start
    idx_start = access_io[0::2] - 1
    idx_end = access_io[1::2] - 1

    return idx_start, idx_end


def teme2AzEl(r, v, time, station_lat, station_lon, station_height):

    #transform r and v coordinates into astropy classes
    teme_p = CartesianRepresentation(r.T*u.km)
    teme_v = CartesianDifferential(v.T*u.km/u.s)
    teme = TEME(teme_p.with_differentials(teme_v), obstime=time)

    #montegancedo coordinates
    station = EarthLocation.from_geodetic(station_lon * u.deg, station_lat * u.deg, station_height * u.m)

    result = teme.transform_to(AltAz(obstime=time, location=station))

    return result.az.value, result.alt.value


def propagate(satellite, starting_date, propagation_time, dT):
    """
    This function propagates satellite using SPG4

    Parameters
    ----------
    satellite: spg4 satellite class
        Satellite orbital information
    starting_date: astropy time class
        Start time for propagation
    propagation_time: float
        Propagation period
    dT: float
        Time step

    Returns
    -------
    teme_p :  nd-array
        Position of satellite in TEME frame
    teme_v :  nd-array
        Velocity of satellite in TEME frame
    time :  1D-array
        Discretization time
    """

    t1 = starting_date.jd
    t2 = t1 + propagation_time

    time_steps = np.linspace(t1, t2, int((t2-t1) * 86400/dT))

    jd = time_steps.astype(int)
    fr = time_steps - time_steps.astype(int)

    error_code, teme_p, teme_v = satellite.sgp4_array(jd, fr)  # in km and km/s

    time = Time(time_steps, format='jd')

    return teme_p, teme_v, time


# def log_TLE(UTC_time):
#     """
#     Open the TLE log and write the TLE of the day

#     Parameters
#     ----------
#     date :  str
#         UTC Date in ISO format: "2021-01-01 12:00:00.000"
#     """

#     #keep only the year-month-day part
#     date = UTC_time[0:10]

#     #Satellite data
#     upmsat_id = '46276'
#     upmsat_TLE = fetch_tle_from_celestrak(upmsat_id)

#     TLE_1 = upmsat_TLE[1]
#     TLE_2 = upmsat_TLE[2]

#     #open log
#     file = open(rel_path + '/logs/TLE_log.txt' , 'a')

#     # write date and TLE
#     file.write(date + '\n')
#     file.write(TLE_1 + '\n')
#     file.write(TLE_2 + '\n')
#     file.close()

#     return 0