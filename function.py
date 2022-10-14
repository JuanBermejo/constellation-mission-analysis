import numpy as np

def Rs(i,e,w):
    """
    Calculates the mean radius under the satellite, Rs [deg]

    Parameters
    ----------
    i : float
        inclination of the orbit [deg]
    e : float
        excentricity of the orbit
    w : float
        argument of the perigee of the orbit [deg]

    Returns
    -------
    Rs : float
        mean radius under the satellite [deg]
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
    Calculates the semimayor axis of an orbit, a [km]

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
