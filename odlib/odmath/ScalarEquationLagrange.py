import numpy as np


def SEL(tau, earth_sun, rhohat, ds):
    # ds = D0, D21, D22, D23
    roots, rhos = np.array([0., 0., 0.]), np.array([0., 0., 0.])
    _tau = tau[2] - tau[0]
    A1 = tau[2] / _tau
    B1 = A1 * (_tau ** 2 - tau[2] ** 2) / 6
    A3 = -tau[0] / _tau
    B3 = A3 * (_tau ** 2 - tau[0] ** 2) / 6
    A = (A1 * ds[1] - ds[2] + A3 * ds[3]) / (-ds[0])
    B = (B1 * ds[1] + B3 * ds[3]) / (-ds[0])
    E = -2 * (np.dot(earth_sun, rhohat))
    F = np.linalg.norm(earth_sun) ** 2
    a = -(A ** 2 + A * E + F)
    b = -(2 * A * B + B * E)
    c = -B ** 2
    coef = np.array([1, 0, a, 0, 0, b, 0, 0, c])
    roots = np.polynomial.polynomial.polyroots(coef)
    roots = roots[np.isreal(np.array(roots))]
    roots = roots[roots > 0]
    print("\nRoots:", roots)
    # Code not working there might be a conversion error


def main():
    tau = [0.15481889055, 0.15481889055, 0.3096377811]
    rhohat2 = [0.8518563498182248, -0.2484702599212149, 0.4610892421311239]
    ds = [-0.0010461861084885213, -0.17297581974209159, -0.17201260125558127, -0.16712421570714076]
    earth_sun = [0.2398478458274071, 0.9065739917845802, 0.392962374977095]
    SEL(tau, earth_sun, rhohat2, ds)


if __name__ == '__main__':
    main()