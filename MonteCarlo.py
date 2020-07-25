from astropy.io import fits
import LueyOD
import numpy as np
import odlib as od
from multiprocessing import Process
import ephem


def read_corr():
    # inputs = ["io/Images/07-03/corr_set1.fits", "io/Images/07-03/corr_set2.fits", "io/Images/07-03/corr_set3.fits",
    #           "io/Images/07-22/corr_set1.fits", "io/Images/07-22/corr_set2.fits",
    #           "io/Images/07-12/corr_set1.fits", "io/Images/07-12/corr_set2.fits",
    #           "io/Images/06-25/corr_set3.fits"]
    inputs = [
              "io/Images/07-03/corr_set1.fits", "io/Images/07-03/corr_set2.fits", "io/Images/07-03/corr_set3.fits"]
    uncert_ra, uncert_dec = np.empty((3)), np.empty((3))

    for c, name in enumerate(inputs):
        table = fits.open(name)[1].data
        # square difference add then sqrt
        mean_dec, mean_ra = 0.0, 0.0
        for i in range(len(table)):
            mean_dec+=(table.index_dec[i] - table.field_dec[i]) ** 2
            mean_ra+=(table.index_ra[i] - table.field_ra[i]) ** 2
        uncert_ra[c],uncert_dec[c] = (mean_ra/len(table)) ** (1/2), (mean_dec/len(table)) ** (1/2)
    for i in range(3): print(od.convert_degrees_HMS(uncert_ra[i]))
    for i in range(3): print(od.convert_degrees_DMS(uncert_dec[i]))


def run(thread):
    _file_name = "io/Observation/1988DZ4.txt"
    elements, iterations = np.array([]), 25000
    times, RA, dec, R, rhohat = LueyOD.parse(_file_name)
    times, RA, dec, R = np.take(times, [0,2,7]), np.take(RA, [0,2,7]), np.take(dec, [0,2,7]), np.take(R, [0,2,7])
    _ra = np.array([np.random.normal(RA[0], 2.166E-4, iterations),
                    np.random.normal(RA[1], 1.3964E-4, iterations),
                    np.random.normal(RA[2], 8.8231E-6, iterations)])
    _dec = np.array([np.random.normal(dec[0], 1.645E-6, iterations),
                     np.random.normal(dec[1], 1.9143E-4, iterations),
                     np.random.normal(dec[2], 2.6904E-4, iterations)])
    for iteration in range(iterations):
        rhohat = np.empty((len(times),), dtype=object)
        for i in range(len(times)): rhohat[i] = np.array(od.convert_RA_dec_rho(_ra[i][iteration], _dec[i][iteration]))
        li = LueyOD.orbit_determination(rhohat, R, times, [0, 1, 2], _input=0)
        if li is None: continue
        li = np.squeeze(li)
        if elements.size == 0: elements = li
        else: elements = np.vstack((elements,li))
        if iteration % 1000 == 0: print(thread, iteration)
    with open('io/MonteCarlo/monte.txt', 'ab') as f: np.savetxt(f, elements, fmt='%.8f')

    # Uncertainty
    # 2.166E-4, 1.3964E-4, 8.8231E-6
    # 1.645E-6, 1.9143E-4, 2.6904E-4


def thread():
    f = open("io/MonteCarlo/monte.txt", 'w')
    f.close()
    print("Running 1,000,000 simulations")
    p, p2, p3, p4 = Process(target=run, args=("CPU Core 1: ", )), Process(target=run, args=("CPU Core 2: ", )), Process(target=run, args=("CPU Core 3: ", )), Process(target=run, args=("CPU Core 4: ", ))
    p.start()
    p2.start()
    p3.start()
    p4.start()
    p.join()
    p2.join()
    p3.join()
    p4.join()
    elements = np.loadtxt("io/MonteCarlo/monte.txt")
    elements = np.reshape(elements, (-1, 10))
    print(np.average(elements, axis=0))
    return elements


def consistency_check():
    # Using pyephem
    obs = ephem.Observer()
    obs.lon, obs.lat, obs.elev = -120.54, 47.002, 490.
    obs.date = "2020/07/22 20:24:40"
    mean_anomaly = od.Ephemeris.calculate_M(1.9057, 3.336483, 2459052.85046, 2459033.75993)
    data = "{name},e,{inclination},{ascending_node},{argument_perihelion},{semimajor_axis},,{eccentricity},{mean_anomaly},{month}/{decimal_days}/{year},2000,,".format(
        name="1988DZ4",
        inclination=19.1805,
        ascending_node=176.1176,
        argument_perihelion=105.3783,
        semimajor_axis=1.9057,
        eccentricity=0.1308,
        mean_anomaly=mean_anomaly,
        month=7,
        decimal_days=22.85046,
        year=2020
    )
    obj = ephem.readdb(data)
    obj.compute(obs)
    ra = od.convert_HMS_degrees(*list(map(float, str(obj.a_ra).split(":"))))
    dec = od.convert_DMS_degrees(*list(map(float, str(obj.a_dec).split(":"))))
    print(ra, dec)


if __name__ == '__main__':
    thread()
