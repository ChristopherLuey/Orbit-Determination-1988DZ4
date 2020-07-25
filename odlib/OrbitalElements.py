import numpy as np
from math import degrees, acos, cos, sin, sqrt
from odlib.Constants import *
from odlib.Conversion import atan2


def calculate_angular_momentum(position, velocity):
    return np.cross(position, velocity)


def calculate_semimajor_axis(position, v):
    # Position and velocity
    return 1/((2/np.linalg.norm(position)) - np.dot(v,v))


def calculate_eccentricity(h, a):
    # Momentum and semimajor axis
    return (1-(np.linalg.norm(h)**2)/a) ** (1/2)


def calculate_orbit_inclination(h):
    # Momentum
    return degrees(acos(h[2]/np.linalg.norm(h)))


def calculate_longitude_of_ascending_node(h, i):
    # Momentum and inclination
    return degrees(atan2((h[0] / (np.linalg.norm(h) * sin(radians(i)))), (-h[1] / (np.linalg.norm(h) * sin(radians(i))))))


def calculate_perihelion(r, i, a, e, r_dot, h, omega):
    # Position, inclination, semimajor axis, eccentricity, velocity, momentum, longitude of ascending node
    v_sin, v_cos = (a * (1 - e ** 2) * np.dot(r, r_dot)) / (np.linalg.norm(h) * np.linalg.norm(r) * e), ((a * (1 - e ** 2) / np.linalg.norm(r)) - 1) / e
    v = atan2(v_sin, v_cos)
    U_cos, U_sin = (r[0] * cos(radians(omega)) + r[1] * sin(radians(omega))) / np.linalg.norm(r), r[2] / (np.linalg.norm(r) * sin(radians(i)))
    U = atan2(U_sin, U_cos)
    return degrees(U) - degrees(v) + 360


def calculate_mean_anomaly(e, E):
    # Eccentricity, eccentric anomaly
    return degrees(E - e * (sin(E)))


def calculate_eccentric_anomaly(e, r, a):
    # eccentricity, position, semimajor axis
    return (acos((1 - np.linalg.norm(r) / a) / e))


def calculate_last_perihelion_passage(M, a, t):
    # mean anomaly, semimajor axis, julian date
    #return -(M*a**(3/2))/k + t
    return -(radians(M) * sqrt(a ** 3)) / k + t


def calculate_n(mu=None, a=None, h=None, e=None, r=None, r_dot=None):
    if mu is not None and a is not None:
        return (mu/a**3) ** (1/2)
    elif h is not None and e is not None and a is not None:
        return h/((a**2) * (1-e**2) ** (1/2))
    elif r is not None and r_dot is not None:
        a,h = calculate_semimajor_axis(r, r_dot), np.linalg.norm(calculate_angular_momentum(r, r_dot))
        return h / ((a ** 2) * (1 - calculate_eccentricity(h, a) ** 2) ** (1 / 2))
    else:
        return 0

