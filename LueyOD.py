from math import radians, degrees

import numpy as np
import zipfile
import os, shutil
if not os.path.isfile('odlib/__init__.py'):
    with zipfile.ZipFile("odlib.zip","r") as zip_ref: zip_ref.extractall("")
    print("Extracting required files...")
    os.remove('odlib.zip')
if not os.path.isfile('io/Observation/LueyInput.txt'):
    with zipfile.ZipFile("io.zip", "r") as zip_ref: zip_ref.extractall("")
    print("Extracting required files...")
    os.remove('io.zip')
try: shutil.rmtree("__MACOSX")
except: pass
import odlib as od
from itertools import combinations


def main():
    _file_name = input("Enter the input file for orbit determination: ")
    #_file_name = "io/Observation/1988DZ4.txt"
    print(_file_name.replace("Observation", "OrbitalElements"))
    output, _file = open(_file_name.replace("Observation", "OrbitalElements"), "w"), open(_file_name, "r")
    pprint("Gathering data from {}\n".format(_file_name), output)
    f = _file.read()

    # Parse
    times, RA, dec, R, rhohat = parse(_file_name)
    pprint(f, file=output, end=""), pprint("\nFound {} recorded data points.".format(len(times)), output)

    # Ask user for which rows to run
    data, text = list(range(len(times))), []
    if len(times) > 3:
        data = input("\nWhich data points would you like to run? Please separate rows using commas.\nNOTE: Rows start at index 0\nNOTE: Enter 'all' for all: ")
        try:
            if data == "all":data = list(range(len(times)))
            else: data=[int(x) for x in list(map(str.strip, data.split(",")))]
        except: pprint("\nWhich data points would you like to run? Please separate rows using commas.\nNOTE: Rows start at index 0\nNOTE: Enter 'all' for all: " + str(data), output), quit()
        pprint("\nWhich data points would you like to run? Please separate rows using commas.\nNOTE: Rows start at index 0\nNOTE: Enter 'all' for all: " + str(data), output)
    data.sort()

    # Run the OD with user inputs
    elements = orbit_determination(rhohat, R, times, data, output, text, _input=True)
    #if elements[0] is not None: print_orbital_elements(elements, text)
    #pprint(text, output)
    output.close()


def orbit_determination(rhohat, R, times, data, output=None, text=None, _input=True):
    elements, combo = [], 0
    for i in list(combinations(np.arange(0, len(times)), 3)):
        exec("i=list(i)\ni.sort()")
        if not all(elem in data for elem in i):
            if text is not None:text = []
            continue
        _times = np.take(times, i)
        _rhohat, _R, tau = np.take(rhohat, i), np.take(R, i), np.array([od.convert_day_to_gaussian(_times[0] - _times[1]), od.convert_day_to_gaussian(_times[2] - _times[1]),od.convert_day_to_gaussian(_times[2] - _times[0])])
        D0, D = od.GaussMethod.calculate_D(_rhohat, _R)
        r2, rho2, roots = od.ScalarEquationLagrange.ScalarEquationLagrange(tau, _R[1], _rhohat[1], np.array([D0, D[1][0], D[1][1], D[1][2]]))
        if _input is True:
            pprint("\nSolving using rows {}, {}, {}:\n\n\tFiltered roots using Scalar Equation of Lagrange:".format(*i), output)
            for c, _ in enumerate(r2): pprint("\t\t({}) r2 = {} AU\t\trho2 = {} AU".format(c + 1, r2[c], rho2[c]), output)
            try: choice = int(input('\n\tEnter manually which root you wish to solve: '));
            except:quit()
        else:
            choice=_input
            if text is not None:
                text.extend(["\nSolving using rows {}, {}, {}:".format(*i), "\n\tFiltered roots using Scalar Equation of Lagrange:", ])
                for c, _ in enumerate(r2): text.append("\t\t({}) r2 = {} AU\t\trho2 = {} AU".format(c + 1, r2[c], rho2[c]))

        # Chooses to test all possible root values
        if choice != 0:
            try: r2 = [r2[choice-1]]
            except:
                if text is not None:
                    text.pop(len(text)-1)
                    text.append("\nInvalid choice");exec("pprint(text,output)\nquit()")
        for r2s in r2:
            #_elements = []
            if text is not None:
                text.append("\n\t#{}; Method of Gauss for root {}:".format(combo, np.where(r2 == r2s)[0] + 1))
            _r2, r2_dot, rho2, text = od.GaussMethod.GaussMethod(_R, _rhohat, r2s, tau, _times, text, D0, D)
            combo+=1
            if _r2 is None or r2_dot is None or rho2 is None or np.isnan(np.sum(_r2)) or np.isnan(np.sum(r2_dot)) or rho2 == float('nan'):
                if text is not None:
                    text.append("\t\tDiverges")
                return None
                #elements.append([0,0,0,0,0,0,0,0,0,0]);continue
            _r2_ecliptic, r2_dot_ecliptic = display_r2(_r2, r2_dot, rho2, text)
            display_orbital_elements(_r2_ecliptic, r2_dot_ecliptic, _times, text, elements)
            #elements.append(_elements)
        if text is not None:
            exec("pprint(text,output)\ntext=[]")
    elements = np.array(elements)
    return elements


def display_r2(_r2, r2_dot, rho2, text):
    _r2_ecliptic, r2_dot_ecliptic = od.convert_cartesian_equatorial_cartesian_ecliptic(_r2), od.convert_cartesian_equatorial_cartesian_ecliptic(r2_dot)
    if text is not None:
        text.extend(["\n\t\tCartesian Equatorial Coordinates:", "\t\t\tr2 = " + str(_r2) + " = " + str(np.linalg.norm(_r2)) + " AU",
                 "\t\t\tr2_dot = " + str(r2_dot) + " = " + str(np.linalg.norm(r2_dot)) + " AU/day", "\t\tCartesian Ecliptic Coordinates:",
                 "\t\t\tr2 = " + str(_r2_ecliptic) + " = " + str(np.linalg.norm(_r2_ecliptic)) + " AU",
                 "\t\t\tr2_dot = " + str(r2_dot_ecliptic) + " = " + str(np.linalg.norm(r2_dot_ecliptic)) + " AU/day",
                 "\n\t\tρ2 =  " + str(rho2) + " AU"])
    return _r2_ecliptic, r2_dot_ecliptic


def display_orbital_elements(_r2_ecliptic, r2_dot_ecliptic, _times, text, elements):
    asteroid = od.Asteroid.Asteroid(_r2_ecliptic, r2_dot_ecliptic, _times[1])
    print_orbital_elements(elements, text, asteroid=asteroid, orbital_elements=asteroid.get_orbital_elements(), _times=_times)


def print_orbital_elements(elements, text=None, asteroid=None, orbital_elements=None, _times=None):
    if asteroid is not None:
        semimajor_axis, eccentricity, inclination, ascending_node, argument_perihelion, mean_anomaly = orbital_elements[0], orbital_elements[1], orbital_elements[2], orbital_elements[3], orbital_elements[4], orbital_elements[5][_times[1]]
        if text is not None:
            text.extend([
                "\n\t\tOrbital Elements:",
                "\t\t\ta = {} AU".format(semimajor_axis),
                "\t\t\te = {}".format(eccentricity),
                "\t\t\ti = {} deg".format(inclination),
                "\t\t\tω = {} deg".format(argument_perihelion),
                "\t\t\tΩ = {} deg".format(ascending_node),
                "\t\t\tE = {} deg".format(degrees(asteroid.get_element("E"))),
                "\t\t\tM = {} deg".format(mean_anomaly),
                "\t\t\tn = {} rad/day".format(asteroid.get_element("n")),
                "\t\t\tLast Perihelion Passage = {} JD".format(asteroid.get_element("perihelion_passage")),
                "\t\t\tP = {} yrs".format(asteroid.get_element("P"))])

        elements.append([semimajor_axis, eccentricity, inclination, argument_perihelion, ascending_node, degrees(asteroid.get_element("E")), mean_anomaly, asteroid.get_element("n"), asteroid.get_element("perihelion_passage"), asteroid.get_element("P")])
    else:
        text.extend([
            "\nAverage Orbital Elements:",
            "a = {} AU".format(elements[0]),
            "e = {}".format(elements[1]),
            "i = {} deg".format(elements[2]),
            "ω = {} deg".format(elements[3]),
            "Ω = {} deg".format(elements[4]),
            "E = {} deg".format(elements[5]),
            "M = {} deg".format(elements[6]),
            "n = {} rad/day".format(elements[7]),
            "Last Perihelion Passage = {} JD".format(elements[8]),
            "P = {} yrs".format(elements[9])
        ])


def HMS(u): return od.convert_HMS_degrees(*list(map(float, u.decode('UTF-8').split(":")))) if u.decode('UTF-8').find(':') != -1 else u


def DMS(u): return od.convert_DMS_degrees(*list(map(float, u.decode('UTF-8').split(":")))) if u.decode('UTF-8').find(':') != -1 else u


def convert(RA, dec): return np.array(od.convert_RA_dec_rho(RA, dec))


def parse(_file_name):
    times, RA, dec, R1, R2, R3 = np.loadtxt(_file_name, unpack=True, converters={1: HMS, 2: DMS},
                                            dtype={'names': ('times', 'RA', 'dec', 'R1', 'R2', 'R3'),
                                                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
    rhohat, R, _R, _rhohat = np.empty((len(times),), dtype=object), np.empty((len(times),), dtype=object), np.array(
        [R1, R2, R3], dtype=object).T, list(
        np.squeeze(np.array([convert(RA[i], dec[i]) for i in range(len(RA))], dtype=float)))
    for i in range(len(times)): rhohat[i], R[i] = np.array(od.convert_RA_dec_rho(RA[i], dec[i])), np.array(
        [_R[i][0], _R[i][1], _R[i][2]])

    return times, RA, dec, R, rhohat


def pprint(message, file, end="\n", sep="\n"):
    (print(*message, file=file, end=end, sep=sep),print(*message, end=end, sep=sep)) if type(message) == list else (print(message, file=file, end=end),print(message, end=end))


if __name__ == '__main__':
    main()
