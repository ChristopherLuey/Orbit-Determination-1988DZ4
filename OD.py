import numpy as np
import odlib as od
# Use builtin datetime library to manage dates
from datetime import datetime


def main():
    # Parse file & gather inputs
    # Analyzing and predicting for asteroid 2002QF15
    _file_name = "input/OrbitalElements/2002QF15.txt"
    _file_name_prediction = "input/Ephemeris/2002QF15.txt"

    f = open(_file_name, "r").readlines()
    solve_date, name, pos, vel, earth_sun = datetime.strptime(f[0].strip(), "%Y/%m/%d %H:%M:%S"), f[1].strip(), np.array(list(map(str.strip, f[2].split("\t")))).astype(np.float), np.array(list(map(str.strip, f[3].split("\t")))).astype(np.float), np.array(list(map(str.strip, f[5].split("\t")))).astype(np.float)

    f = open(_file_name_prediction, "r").readlines()
    predictions = []
    for i in range(0,len(f), 2): predictions.append((datetime.strptime(f[i].strip(), "%Y/%m/%d %H:%M:%S"), np.array(list(map(str.strip, f[i+1].split("\t")))).astype(np.float)))

    # Create an asteroid object
    asteroid = od.Asteroid.Asteroid(pos, vel, solve_date)
    asteroid.set_name(name)

    # Calculate the RA/Dec for all desired prediction times
    for prediction in predictions:
        prediction_datetime, earth_sun_ = prediction[0], prediction[1]
        asteroid.set_earth_sun_vector(prediction_datetime, earth_sun_)
        prediction_date = od.calculate_julian_date(datetime=prediction_datetime)
        print("\nERROR\n\tDetected incorrect M...\n\tOverriding with expected M")
        asteroid.set_mean_anomaly(prediction[0], M=158.5663308495127)
        ra, dec = asteroid.calculate_RA_dec(prediction_datetime)
        _ra, _dec = od.Ephemeris.calculate_RA_dec_ephem(prediction_datetime, asteroid)
        print("\nRA/Dec Prediction:\n\tJulian Date:", prediction_date, "\n\tDate (mm:dd:yyyy HH:MM:SS):", prediction_datetime.strftime("%m/%d/%Y %H:%M:%S"), "\n\tEarth Sun Vector:", earth_sun_, "\n\tCalculated:", ra, dec, "\n\tPyephem Expected:", _ra, _dec)

    # Gather orbital elements from asteroid object
    orbital_elements = asteroid.get_orbital_elements()

    print("\n\nThe following orbital elements were calculated using:\n\tJulian Date:", asteroid.get_element("calculated_date"), "\n\tDate (mm/dd/yyyy hh:mm:ss):", solve_date.strftime("%m/%d/%Y %H:%M:%S"))
    semimajor_axis, eccentricity, inclination, ascending_node, argument_perihelion, mean_anomaly = orbital_elements[0], orbital_elements[1], orbital_elements[2], orbital_elements[3], orbital_elements[4], orbital_elements[5][solve_date]
    perihelion_passage = asteroid.get_element("perihelion_passage")
    semimajor_axis_expected, eccentricity_expected, inclination_expected, ascending_node_expected, argument_perihelion_expected, mean_anomaly_expected, perihelion_passage_expected = 1.056679892399276, 0.3442561302111088, 25.15346367650418, 236.2549873540567, 255.5046578188225, 139.5121254404565, 2458158.720843222924
    print("\nSemimajor Axis:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(semimajor_axis, semimajor_axis_expected, 100*abs((semimajor_axis - semimajor_axis_expected)/semimajor_axis_expected)))
    print("\nEccentricity:\n\tCalculated Value: {}\n\tExpected Value:{}\n\tError: {}%".format(eccentricity, eccentricity_expected, 100*abs((eccentricity - eccentricity_expected)/eccentricity_expected)))
    print("\nInclination:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(inclination, inclination_expected, 100*abs((inclination - inclination_expected)/inclination_expected)))
    print("\nAscending Node:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(ascending_node, ascending_node_expected, 100*abs((ascending_node - ascending_node_expected)/ascending_node_expected)))
    print("\nArgument of Perihelion:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(argument_perihelion, argument_perihelion_expected, 100*abs((argument_perihelion - argument_perihelion_expected)/argument_perihelion_expected)))
    print("\nMean Anomaly:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(mean_anomaly, mean_anomaly_expected, 100*abs((mean_anomaly - mean_anomaly_expected)/mean_anomaly_expected)))
    print("\nLast Perihelion Passage:\n\tCalculated Value: {}\n\tExpected Value: {}\n\tError: {}%".format(perihelion_passage, perihelion_passage_expected, 100*abs((perihelion_passage - perihelion_passage_expected)/perihelion_passage_expected)))


if __name__ == '__main__':

    # RA/Dec calculation doesn't work. I know that my calculation for mean anomaly on Aug 3 is incorrect, not sure if
    # it's a units problem. When applying the correct mean anomaly, the RA and dec I calculate is incorrect. I
    # believe that there is an issue with radians and degrees since my E value is correct
    main()
