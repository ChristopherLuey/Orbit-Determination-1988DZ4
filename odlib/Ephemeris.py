import numpy as np
from math import degrees, radians, sin, cos, asin, acos
import ephem
from odlib.Conversion import convert_HMS_degrees, convert_DMS_degrees, calculate_decimal_days
from odlib.Constants import *
from odlib.odmath import NewtonRaphson


def calculate_M(a, M, t, t2):
    # Semimajor axis, mean anomaly t0 (aug 3), t2 (july 13)
    print("\nCalculated M for Aug 3:", (t - t2)*((k**2/a**3)**(1/2)) + M)
    print("Expected M for Aug 3:", "158.5663308495127")
    # Something is going wrong here
    return (t - t2)*((k**2/a**3)**(1/2)) + M


def calculate_E(M, e):
    # M at ephemeris time, eccentricity
    return NewtonRaphson.newton_raphson(M, e, 1e-10)


def calculate_cartesian_coordinates(a, E, e):
    # Semimajor axis, E from Kepler's, eccentricity
    # Calculates position vector in cartesian
    return np.array([[a*cos(E) - a*e],[a*sin(E)*(1-e**2)**(1/2)],[0]])


def convert_cartesian_ecliptic(r, OM, i, omega):
    # position cartesian, longitude of ascending node, inclination, argument of periapsis
    OM, i, omega = radians(OM), radians(i), radians(omega)
    rot_a = np.array([[cos(OM), -sin(OM), 0], [sin(OM), cos(OM), 0], [0, 0, 1]])
    rot_b = np.array([[1,0, 0], [0, cos(i), -sin(i)], [0, sin(i), cos(i)]])
    rot_c = np.array(([cos(omega), -sin(omega), 0], [sin(omega), cos(omega), 0], [0, 0, 1]))
    return np.dot(np.dot(rot_a, rot_b), np.dot(rot_c, r))


def convert_ecliptic_cartesian(r):
    # position ecliptic
    rot_a = np.array([[1, 0, 0], [0, (cos(radians(EPSILON))), -(sin(radians(EPSILON)))], [0, (sin(radians(EPSILON))), (cos(radians(EPSILON)))]])
    return np.dot(rot_a, r)


def calculate_rho(r, R):
    # position r in equatorial, earth to sun vector
    # Returns rho, rhohat
    rho = np.array([r[0] + R[0], r[1] + R[1], r[2] + R[2]])
    return rho, rho/np.linalg.norm(rho)


def calculate_RA_dec(rhohat):
    dec = asin(rhohat[2])
    ra = acos(rhohat[0] / cos(dec))
    return degrees(ra), degrees(dec)


def calculate_RA_dec_ephem(prediction_date, asteroid):
    # Using pyephem
    obs = ephem.Observer()
    obs.lon, obs.lat, obs.elev = "0:0:0", "0:0:0", 0.
    obs.date = prediction_date.strftime("%Y/%m/%d %H:%M:%S")
    data = "{name},e,{inclination},{ascending_node},{argument_perihelion},{semimajor_axis},,{eccentricity},{mean_anomaly},{month}/{decimal_days}/{year},2000,,".format(
        name=asteroid.get_name(),
        inclination=asteroid.get_element("inclination"),
        ascending_node=asteroid.get_element("ascending_node"),
        argument_perihelion=asteroid.get_element("argument_perihelion"),
        semimajor_axis=asteroid.get_element("semimajor_axis"),
        eccentricity=asteroid.get_element("eccentricity"),
        mean_anomaly=asteroid.get_element("mean_anomaly")[prediction_date],
        month=prediction_date.month,
        decimal_days=calculate_decimal_days(datetime=prediction_date),
        year=prediction_date.year
    )
    obj = ephem.readdb(data)
    obj.compute(obs)
    ra = convert_HMS_degrees(*list(map(float, str(obj.a_ra).split(":"))))
    dec = convert_DMS_degrees(*list(map(float, str(obj.a_dec).split(":"))))
    asteroid.set_RA(prediction_date, ra)
    asteroid.set_dec(prediction_date, dec)
    return ra, dec
