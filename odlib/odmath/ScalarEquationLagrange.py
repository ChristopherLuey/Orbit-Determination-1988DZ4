import numpy as np


def ScalarEquationLagrange(tau, earth_sun2, rhohat2, ds):
    # ds = D0, D21, D22, D23
    # tau = tau1, tau3, tau
    # roots = r2
    roots, rho2 = np.array([]), np.array([])
    A1, A3 = tau[1] / tau[2], -tau[0] / tau[2]
    B1, B3 = A1 * (tau[2] ** 2 - tau[1] ** 2) / 6, A3 * (tau[2] ** 2 - tau[0] ** 2) / 6
    A, B = (A1 * ds[1] - ds[2] + A3 * ds[3]) / (-ds[0]), (B1 * ds[1] + B3 * ds[3]) / (-ds[0])
    E, F = -2 * (np.dot(earth_sun2, rhohat2)), np.linalg.norm(earth_sun2) ** 2
    a, b, c = -(A ** 2 + A * E + F), -(2 * A * B + B * E), -B ** 2
    roots = np.roots(np.array([1, 0, a, 0, 0, b, 0, 0, c]))
    roots = roots[np.isreal(roots)][roots[np.isreal(roots)] > 0]
    for i in roots: rho2 = np.append(rho2, A + B/(i**3))
    return roots, rho2


def test():
    tau2 = [-0.15481889055, 0.15481889055, 0.3096377811]
    earth_sun2 = [-0.2398478458274071, 0.9065739917845802, 0.392962374977095]
    rhohat2 = [-0.8518563498182248, -0.2484702599212149, 0.4610892421311239]
    ds = [-0.0010461861084885213, -0.17297581974209159, -0.17201260125558127, -0.16712421570714076]
    ScalarEquationLagrange(tau2, earth_sun2, rhohat2, ds)


if __name__ == '__main__':
    test()