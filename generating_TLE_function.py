import numpy as np

def Rs(i,e,w):
    """
    Calculates the mean radius under the satellite, Rs [km]

    Parameters
    ----------
    i : float
        inclination of the orbit [deg]
    e : float
        eccentricity of the orbit [-]
    w : float
        argument of the perigee of the orbit [deg]

    Returns
    -------
    Rs : float
        mean radius under the satellite [km]
    """

    coef_matrix = np.matrix(
        [[0, 0, 0,  6.377788E6],
        [1, 1, 2,  -9.972427E-2],
        [0, 0, 1,  2.836106E1],
        [1, 1, 3,  7.383957E-4],
        [0, 0, 2,  -4.915346E-1],
        [1, 2, 0,  9.150674E1],
        [0, 0, 3,  -4.017214E-3],
        [1, 2, 1,  -1.33747E0],
        [0, 0, 4,  1.273759E-4],
        [1, 2, 2,  1.155062E-3],
        [0, 0, 5,  -5.643483E-7],
        [1, 3, 0,  4.645793E0],
        [0, 1, 0,  2.120874E3],
        [1, 3, 1,  -1.909859E0],
        [0, 1, 1,  -1.886565E2],
        [1, 4, 0,  5.290634E1],
        [0, 1, 2,  5.241726E0],
        [2, 0, 0,  -3.926748E0],
        [0, 1, 3,  -3.84173E-2],
        [2, 0, 1,  1.703196E-2],
        [0, 1, 4,  -2.268633E-6],
        [2, 0, 2,  -8.207037E-6],
        [0, 2, 0,  -3.447518E3],
        [2, 0, 3,  2.674298E-8],
        [0, 2, 1,  8.553227E1],
        [2, 1, 0,  4.368207E0],
        [0, 2, 2,  -4.719285E0],
        [2, 1, 1,  -9.744218E-2],
        [0, 2, 3,  3.438466E-2],
        [2, 1, 2,  -8.298287E-6],
        [0, 3, 0,  9.952275E3],
        [2, 2, 0,  2.116514E-1],
        [0, 3, 1,  1.634122E2],
        [2, 2, 1,  -4.547244E-4],
        [0, 3, 2,  5.259614E-2],
        [2, 3, 0,  -9.768337E-2],
        [0, 4, 0,  -1.695787E4],
        [3, 0, 0,  1.323687E-3],
        [0, 4, 1,  -8.008678E1],
        [3, 0, 1,  -1.327782E-4],
        [0, 5, 0,  8.501297E3],
        [3, 0, 2,  5.032226E-8],
        [1, 0, 0,  2.948006E1],
        [3, 1, 0,  -3.319618E-2],
        [1, 0, 1,  -1.151884E0],
        [3, 1, 1,  7.264131E-4],
        [1, 0, 2,  1.653916E-2],
        [3, 2, 0,  -5.846808E-4],
        [1, 0, 3,  -1.276845E-4],
        [4, 0, 0,  4.655915E-4],
        [1, 0, 4,  3.521573E-8],
        [4, 0, 1,  6.438237E-8],
        [1, 1, 0,  -1.676297E2],
        [4, 1, 0,  4.999741E-6],
        [1, 1, 1,  6.645437E0],
        [5, 0, 0,  -2.08623E-6]]
    )

    rows, colums = coef_matrix.shape

    Rs = 0
    for j in range(rows):
        row_rs = coef_matrix[j,3]*i**coef_matrix[j,0]*e**coef_matrix[j,1]*w**coef_matrix[j,2]
        Rs += row_rs

    return Rs/1000

def semimayor_axis(Ap, e, Rs):
    """
    Calculates the semi-major axis of an orbit, a [km]

    Parameters
    ----------
    Ap : float
        Altitud of the apogee [km]
    e : float
        eccentricity of the orbit [-]
    Rs : float
        mean radius under the satellite [km]

    Returns
    -------
    a : float
        semi-major axis of the orbit [km]
    """

    a = (Ap + Rs)/(1-e)
    return a

def MeanMotion(a):
    """
    Calculates the mean motion of a satellite, MM [rev/day]

    Parameters
    ----------
    a : float
        semi-major axis of the orbit [km]

    Returns
    -------
    MM : float
        mean motion of the satellite [rev/day]
    """

    T_day = 86400 # solar day [s]
    mu = 3.986E5 # standard gravitational parameter of Earth [km^3/s^2]
    MM = T_day/(2*np.pi)*np.sqrt(mu/a**3)
    return MM

def checksum_mod10(variable):
    counter = 0
    for c in variable:  
        if c == "-":
            counter += 1
        elif c.isnumeric():
            counter += int(c)
    checksum = str(counter % 10)
    return checksum

def Generate_TLE(Inclination, RAAN, Eccentricity, AOP, MeanAnomaly, MeanMotion, MeanMotionDot, 
MeanMotionDDot, BSTAR, Epoch, Satellite_catalog_number,
Revolution_number, Classification, International_designator, Element_set_number):
    """
    Generates realistics Two-Line Elements (TLE) for notional systems
    Reference: ROCKWOOD, Troy; STEEGER, Greg; STEIN, Matthew. Generating Realistic Two-Line Elements for Notional 
    Space Vehicles and Constellations. arXiv preprint arXiv:2203.04204, 2022.

    Parameters
    ----------
    Inclination: string
        inclination of the orbit [deg]
    RAAN: string
        right ascension of the ascending node [deg]
    Eccentricity:string
        eccentricity of the orbit [-]
    AOP: string 
        argument of perigee [deg]
    MeanAnomaly: string
        mean anomaly of the orbit [deg]
    MeanMotion: string
        mean motion of the satellite [rev/day]
    MeanMotionDot: string (+.XXXXXXXX)
        first time derivative of the mean motion []
    MeanMotionDDot: string (+XXXXX-X)
        second time derivative of the mean motion []
        leading decimal point assumed
    BSTAR: string (+XXXXX-X)
        drag term []
        leading decimal point assumed
    Epoch: string (XXYYY.YYYYYYYY)
        XX: last two digits of year
        YYY.YYYYYYYY: day of the year and fractional portion of the day
    Satellite_catalog_number: string
        -
    Revolution_number: string
        Revolution number at epoch [revs]
    Classification: string
        U for unclissified
    International_designator: string (XXYYYZZZ)
        XX: last two digits of launch year
        YYY: launch number of the year
        ZZZ: piece of the launch
    Element_set_number: string
        - 
    
    Notes: 
    * The classification, international designator value, element set number and revolution number can
    all be chosen at will, do not affect the simulation
    * The satellite catalog number should not conflict with other satellites simulations
    * The ephemeris type has been set to 0 
    * The epoch can be chosen at will but with some considerations. If working with notional systems, then select similar epochs.
    If mixing notional with real satellites, choose an epoch close to the real systems of interest. 
    * The inclination, RAAN, eccentricity, AOP, and mean anomaly are all chosen according to the desired simulation.
    * The first and second time derivatives of the mean motion and, the drag term can be derived from existing TLE. This terms are 
    important to obtain realistic TLE. 

    Returns
    -------
    TLE : string
        TLE of a notional system
    """

    TLE1 = "1" + " " + Satellite_catalog_number + Classification + " " + International_designator + " " + Epoch + " "
    TLE1 += MeanMotionDot + " " + MeanMotionDDot + " " + BSTAR + " " + "0" + Element_set_number 

    TLE2 = "2" + " " + Satellite_catalog_number + " " + Inclination + " " + RAAN + " " + Eccentricity + " " + AOP + " " 
    TLE2 += MeanAnomaly + " " + MeanMotion + " " + Revolution_number
    
    TLE1 += checksum_mod10(TLE1)
    TLE2 += checksum_mod10(TLE2)

    TLE = (TLE1, TLE2)
    return(TLE)



