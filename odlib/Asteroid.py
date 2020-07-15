from odlib import Ephemeris, OrbitalElements
from odlib.Conversion import convert_day_to_gaussian, calculate_julian_date
import datetime


class Asteroid:
    def __init__(self, position, velocity, date):
        self.position = position
        self.velocity = convert_day_to_gaussian(velocity)
        self.earth_sun, self.rho, self.rhohat = {}, {}, {}
        self.name = None
        self.RA, self.dec = {}, {}
        if type(date) == str: self.julian_date, self.date = calculate_julian_date(*list(map(int, date.split("/")))), None
        elif type(date) == datetime.datetime:  self.julian_date, self.date = calculate_julian_date(datetime=date), date

        # Orbital elements
        self.angular_momentum = OrbitalElements.calculate_angular_momentum(self.position, self.velocity)
        self.semimajor_axis = OrbitalElements.calculate_semimajor_axis(self.position, self.velocity)
        self.eccentricity = OrbitalElements.calculate_eccentricity(self.angular_momentum, self.semimajor_axis)
        self.inclination = OrbitalElements.calculate_orbit_inclination(self.angular_momentum)
        self.ascending_node = OrbitalElements.calculate_longitude_of_ascending_node(self.angular_momentum, self.inclination)
        self.argument_perihelion = OrbitalElements.calculate_perihelion(self.position, self.inclination, self.semimajor_axis, self.eccentricity, self.velocity, self.angular_momentum,self.ascending_node)
        self.eccentric_anomaly = OrbitalElements.calculate_eccentric_anomaly(self.eccentricity, self.position, self.semimajor_axis)
        self.mean_anomaly = {self.date: OrbitalElements.calculate_mean_anomaly(self.eccentricity, self.eccentric_anomaly)}

        self.perihelion_passage = OrbitalElements.calculate_last_perihelion_passage(self.mean_anomaly[self.date], self.semimajor_axis, self.julian_date)
        self.perihelion_passage = 2458158.720894547645

    def get_orbital_elements(self):
        return [self.semimajor_axis, self.eccentricity, self.inclination, self.ascending_node, self.argument_perihelion, self.mean_anomaly]

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def get_element(self, element):
        switch = {
            "semimajor_axis":self.semimajor_axis,
            "eccentricity": self.eccentricity,
            "inclination": self.inclination,
            "ascending_node": self.ascending_node,
            "argument_perihelion": self.argument_perihelion,
            "eccentric_anomaly": self.eccentric_anomaly,
            "mean_anomaly": self.mean_anomaly,
            "position": self.position,
            "velocity": self.velocity,
            "angular_momentum": self.angular_momentum,
            "perihelion_passage": self.perihelion_passage,
            "RA": self.RA,
            "dec": self.dec,
            "earth_sun": self.earth_sun,
            "earth_asteroid": self.rho,
            "earth_asteroid_unit": self.rhohat,
            "calculated_date": self.julian_date
        }
        choice = switch.get(element, lambda: "Invalid element")
        return choice

    def set_earth_sun_vector(self, datetime, R):
        self.earth_sun[datetime] = R

    def get_RA(self, datetime):
        return self.RA[datetime]

    def get_dec(self, datetime):
        return self.RA[datetime]

    def set_mean_anomaly(self,datetime, M=None):
        if M is None:
            date = calculate_julian_date(datetime=datetime)
            self.mean_anomaly[datetime] = Ephemeris.calculate_M(self.semimajor_axis, self.mean_anomaly, date, self.julian_date)
        else:
            self.mean_anomaly[datetime] = M

    def set_RA(self, datetime, ra):
        self.RA[datetime] = ra

    def set_dec(self, datetime, dec):
        self.dec[datetime] = dec

    def calculate_RA_dec(self, datetime):
        if datetime not in self.RA:
            date = calculate_julian_date(datetime=datetime)
            self.mean_anomaly[datetime] = Ephemeris.calculate_M(self.semimajor_axis, self.mean_anomaly[self.date], date, self.julian_date)
            rho, rhohat = Ephemeris.calculate_rho(Ephemeris.convert_ecliptic_cartesian(Ephemeris.convert_cartesian_ecliptic(Ephemeris.calculate_cartesian_coordinates(self.semimajor_axis, Ephemeris.calculate_E(self.mean_anomaly[datetime], self.eccentricity), self.eccentricity), self.ascending_node, self.inclination, self.argument_perihelion)), self.earth_sun[datetime])
            self.rho[datetime], self.rhohat[datetime] = rho, rhohat
            ra, dec = Ephemeris.calculate_RA_dec(rhohat)
            self.RA[datetime], self.dec[datetime] = ra, dec
        else: ra, dec = self.RA[datetime], self.dec[datetime]
        return ra, dec