import numpy as np
import odlib as od
from itertools import combinations


def main():
    #_file_name = input("Enter the input file for orbit determination: ")
    _file_name = "io/Observation/2020test.txt"
    print(_file_name.replace("Observation", "OrbitalElements"))
    output, _file = open(_file_name.replace("Observation", "OrbitalElements"), "w"), open(_file_name, "r")
    pprint("Gathering data from {}\n".format(_file_name), output)

    # Read user inputs define variables
    f = _file.readlines()
    R, times, radec, radec_deg, rhohat = np.empty((len(f),), dtype=object), np.array([]), np.empty((0,2), dtype=float), np.empty((0,2)), np.empty((len(f),), dtype=object)

    # Convert user inputs to proper format
    for _, i in enumerate(f):
        R[_], times, radec = np.array(list(map(lambda x: float(x.strip()), i.split(" ")[3:6]))), np.append(times, float(i.split(" ")[0])), np.vstack((radec, np.array([i.split(" ")[1:3]])))
        pprint(i, file=output, end="")

    pprint("\nFound {} recorded data points.".format(len(times)), output)

    # Convert RA/Dec to degrees
    for c, i in enumerate(radec):
        l = np.array([od.convert_HMS_degrees(*list(map(float, i[0].split(":")))), od.convert_DMS_degrees(*list(map(float, i[1].split(":"))))])
        radec_deg = np.vstack([radec_deg, l])
        rhohat[c] = np.array(od.convert_RA_dec_rho(radec_deg[c][0], radec_deg[c][1]))

    # Ask user for which rows to run
    data, text = list(range(len(times))), []
    if len(times) > 3:
        data = input("\nWhich data points would you like to run? Please separate rows using commas.\nNOTE: Rows start "
                     "at index 0\nNOTE: Enter 'all' for all: ")
        if data == "all": data = list(range(len(times)))
        else: data = [int(x) for x in list(map(str.strip, data.split(",")))]
        pprint("\nWhich data points would you like to run? Please separate rows using commas.\nNOTE: Rows start at "
               "index 0\nNOTE: Enter 'all' for all: " + str(data), output)
    data.sort()

    # Run the OD with user inputs
    elements = orbit_determination(rhohat, R, times, output, data, text, _input=True)
    if elements[0] is not None: print_orbital_elements(elements, text)
    pprint(text, output)
    output.close()


def orbit_determination(rhohat, R, times, output, data, text, _input=True):
    elements = []
    for i in list(combinations(np.arange(0, len(times)), 3)):
        i = list(i)
        i.sort()
        if not all(elem in data for elem in i):
            text = []
            continue
        _times = np.take(times, i)
        _rhohat, _R, tau = np.take(rhohat, i), np.take(R, i), np.array(
            [od.convert_day_to_gaussian(_times[0] - _times[1]), od.convert_day_to_gaussian(_times[2] - _times[1]),
             od.convert_day_to_gaussian(_times[2] - _times[0])])
        D0, D = od.GaussMethod.calculate_D(_rhohat, _R)
        r2, rho2, roots = od.ScalarEquationLagrange.ScalarEquationLagrange(tau, _R[1], _rhohat[1], np.array([D0, D[1][0], D[1][1], D[1][2]]))
        if _input is True:
            pprint("\nSolving using rows {}, {}, {}:\n\n\tFiltered roots using Scalar Equation of Lagrange:".format(*i), output)
            for c, _ in enumerate(r2): pprint("\t\t({}) r2 = {} AU\t\trho2 = {} AU".format(c + 1, r2[c], rho2[c]), output)
            try: choice = int(input("\n\tEnter manually which root you wish to solve: "))
            except: quit()
        else:
            choice=_input
            text.extend(["\nSolving using rows {}, {}, {}:".format(*i), "\n\tFiltered roots using Scalar Equation of Lagrange:", ])
            for c, _ in enumerate(r2): text.append("\t\t({}) r2 = {} AU\t\trho2 = {} AU".format(c + 1, r2[c], rho2[c]))
        text.append("\n\tMethod of Gauss for root {}:".format(choice))

        # Chooses to test all possible root values
        if choice == 0:
            _elements = []
            for r2s in r2:
                _r2, r2_dot, rho2 = od.GaussMethod.GaussMethod(_R, _rhohat, r2s, tau, _times, text, D0, D)
                _r2_ecliptic, r2_dot_ecliptic = display_r2(_r2, r2_dot, rho2, text)
                display_orbital_elements(_r2_ecliptic, r2_dot_ecliptic, _times, text, _elements)
            elements.append(list(np.average(_elements, axis=0)))
        # Chooses specific root value to test
        else:
            try:
                _r2, r2_dot, rho2 = od.GaussMethod.GaussMethod(_R, _rhohat, r2[choice - 1], tau, _times, text, D0, D)
                _r2_ecliptic, r2_dot_ecliptic = display_r2(_r2, r2_dot, rho2, text)
                display_orbital_elements(_r2_ecliptic, r2_dot_ecliptic, _times, text, elements)
            except:
                text.pop(len(text)-1)
                text.append("\nInvalid choice")
                pprint(text, output)
                quit()
        pprint(text, output)
        text = []
    elements = list(np.average(elements, axis=0))
    return elements


def display_r2(_r2, r2_dot, rho2, text):
    _r2_ecliptic, r2_dot_ecliptic = od.convert_cartesian_equatorial_cartesian_ecliptic(_r2), od.convert_cartesian_equatorial_cartesian_ecliptic(r2_dot)
    text.extend(["\n\t\tCartesian Equatorial Coordinates:", "\t\t\tr2 = " + str(_r2) + " = " + str(np.linalg.norm(_r2)) + " AU",
                 "\t\t\tr2_dot = " + str(r2_dot) + " = " + str(np.linalg.norm(r2_dot)) + " AU/day", "\t\tCartesian Ecliptic Coordinates:",
                 "\t\t\tr2 = " + str(_r2_ecliptic) + " = " + str(np.linalg.norm(_r2_ecliptic)) + " AU",
                 "\t\t\tr2_dot = " + str(r2_dot_ecliptic) + " = " + str(np.linalg.norm(r2_dot_ecliptic)) + " AU/day",
                 "\n\t\tρ2 =  " + str(rho2) + " AU"])
    return _r2_ecliptic, r2_dot_ecliptic


def display_orbital_elements(_r2_ecliptic, r2_dot_ecliptic, _times, text, elements):
    asteroid = od.Asteroid.Asteroid(_r2_ecliptic, r2_dot_ecliptic, _times[1])
    orbital_elements = asteroid.get_orbital_elements()
    print_orbital_elements(elements, text, asteroid=asteroid, orbital_elements=orbital_elements, _times=_times)


def print_orbital_elements(elements, text=None, asteroid=None, orbital_elements=None, _times=None):
    if asteroid is not None:
        semimajor_axis, eccentricity, inclination, ascending_node, argument_perihelion, mean_anomaly = orbital_elements[0], orbital_elements[1], orbital_elements[2], orbital_elements[3], orbital_elements[4], orbital_elements[5][_times[1]]
        text.extend([
            "\n\t\tOrbital Elements:",
            "\t\t\ta = {} AU".format(semimajor_axis),
            "\t\t\te = {}".format(eccentricity),
            "\t\t\ti = {} deg".format(inclination),
            "\t\t\tω = {} deg".format(argument_perihelion),
            "\t\t\tΩ = {} deg".format(ascending_node),
            "\t\t\tE = {} deg".format(asteroid.get_element("E")),
            "\t\t\tM = {} deg".format(mean_anomaly),
            "\t\t\tn = {} rad/day".format(asteroid.get_element("n")),
            "\t\t\tLast Perihelion Passage = {} JD".format(asteroid.get_element("perihelion_passage")),
            "\t\t\tP = {} yrs".format(asteroid.get_element("P"))
        ])
        elements.append([semimajor_axis, eccentricity, inclination, argument_perihelion, ascending_node, asteroid.get_element("E"), mean_anomaly, asteroid.get_element("n"), asteroid.get_element("perihelion_passage"), asteroid.get_element("P")])
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


def pprint(message, file, end="\n", sep="\n"):
    if type(message) == list:
        print(*message, file=file, end=end, sep=sep)
        print(*message, end=end, sep=sep)
    else:
        print(message, file=file, end=end)
        print(message, end=end)


if __name__ == '__main__':
    main()
