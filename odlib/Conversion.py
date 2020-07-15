from odlib.Constants import k
from math import acos, asin, pi


def convert_day_to_gaussian(days):
    return days/k


def convert_DMS_degrees(degrees, minutes, seconds):
    # Handle a negative angle
    a=1
    if degrees < 0.0 or degrees == -0.0: a = -1
    # Perform angle conversion
    b = abs(degrees) + (minutes/60) + (seconds/3600)
    # Return result
    return b*a


def convert_HMS_degrees(hours, minutes, seconds):
    # Handle a negative angle
    a=1
    if hours < 0.0 or hours == -0.0: a = -1
    # Perform angle conversion
    b = (abs(hours) + (minutes/60) + (seconds/3600))*15
    # Return result
    return b*a


def convert_degrees_DMS(ang):
    d = int(ang)
    m = 60*(ang - d)
    s, m = (m - int(m)) * 60, int(m)
    return d, m, s


def convert_degrees_HMS(ang):
    h = int(ang/15)
    m = 60*(ang-h*15)/15
    s, m = (m - int(m)) * 60, int(m)
    return h, m, s


def atan2(sine, cosine):
    if cosine > 0 and sine > 0: return asin(sine)
    if cosine < 0 < sine: return acos(cosine)
    if cosine < 0 and sine < 0: return pi - asin(sine)
    if cosine > 0 > sine: return 2 * pi + asin(sine)


def calculate_julian_date(d=None, m=None, y=None, datetime=None):
    if datetime is not None:
        return 367*datetime.year - int(7*(datetime.year + int((datetime.month+9)/12)) / 4) + int(275*datetime.month/9) + calculate_decimal_days(datetime=datetime) + 1721013.5
    else:
        return 367*y - int(7*(y + int((m+9)/12)) / 4) + int(275*m/9) + d + 1721013.5


def calculate_decimal_days(hours=None, minutes=None, seconds=None, days=None, datetime=None):
    if datetime is not None:
        return datetime.second * 1.1574074e-5 + datetime.minute * 0.00069444444 + datetime.hour * 0.041666667 + datetime.day
    else:
        return seconds * 1.1574074e-5 + minutes * 0.00069444444 + hours * 0.041666667 + days

