from generating_TLE_function import *

def generate_FOSSA_TLE(Inclination, RAAN, Eccentricity, AOP, MeanAnomaly, AltitudePerigee, Epoch, Satellite_catalog_number):
    """
    Generates realistics Two-Line Elements (TLE) for notional systems in the FOSSA constellation

    Parameters
    ----------
    Inclination: float
        inclination of the orbit [deg]
    RAAN: float
        right ascension of the ascending node [deg]
    Eccentricit: float
        eccentricity of the orbit [-]
    AOP: float
        argument of perigee [deg]
    MeanAnomaly: float
        mean anomaly of the orbit [deg]
    SemimajorAxis: float
        semi-major axis of the orbit [km]
    Epoch: string (XXYYY.YYYYYYYY)
        XX: last two digits of year
        YYY.YYYYYYYY: day of the year and fractional portion of the day
    Satellite_catalog_number: string
        -

    Notes:
    * The satellite catalog number should not conflict with other satellites simulations
    * The epoch can be chosen at will but with some considerations. If working with notional systems, then select similar epochs.
    If mixing notional with real satellites, choose an epoch close to the real systems of interest.

    Returns
    -------
    FOSSA_TLE : string
        TLE of a notional system for the FOSSA constellation
    """

    # Function: MeanMotion(a) --> for a function semimayor_axis(Ap, e, Rs) --> for Rs function Rs(i, e, AOP)
    MM = MeanMotion(semimayor_axis(AltitudePerigee - Rs(Inclination, Eccentricity, AOP), Eccentricity, Rs(Inclination, Eccentricity, AOP)))

    # The values taken for the first and second time derivative of mean motion and the drag term are the mean values obtained from
    # the evolution of the paremeters from the 1st to the 31th of October 2022 calculated in FOSSA_TLE_evolution
    MeanMotionDot = "+.00030168" # mean value: 0.000301682
    MeanMotionDDot = "00000+0" # mean value: 0.0
    BSTAR = "14049-2" # mean value: 0.0014049229
    MeanMotionDot = 0.00030168
    MeanMotionDDot = 0.
    BSTAR = 0.0014049229

    # The following parameters do not affect the simulation.
    # They have been chosen according to the TLE from FOSSASAT-2E2 (09/11/2022)
    Revolution_number = "4545"
    Classification =  "U"
    International_designator = "22002A"
    Element_set_number = 999

    # The input arguments for the function Generate_TLE are strings
    # FOSSA_TLE = Generate_TLE(str(Inclination), str(RAAN), str(Eccentricity), str(AOP), str(MeanAnomaly), str(MM), MeanMotionDot,
    # MeanMotionDDot, BSTAR, Epoch, Satellite_catalog_number,
    # Revolution_number, Classification, International_designator, Element_set_number)

    FOSSA_TLE = Generate_TLE_2(Inclination, RAAN, Eccentricity, AOP, MeanAnomaly, MM, MeanMotionDot,
    MeanMotionDDot, BSTAR, Epoch, Satellite_catalog_number,
    Revolution_number, Classification, International_designator, Element_set_number)

    return FOSSA_TLE