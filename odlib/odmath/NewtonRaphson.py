from math import sin, cos, inf
import numpy as np
import odlib as od

def _test():
	M, e, tol, dif = 0.42, 0.8, 1e-10, float(inf)
	E, iterations = M, 1
	while dif > tol:
		E = _newton_raphson_mean_anomaly(E, e, M)
		_M = E - e * sin(E)
		dif = abs(M - _M)
		iterations+=1

	print("E =", E)
	print("Convergence Parameter =", tol)
	print("Iterations =", iterations)


def _newton_raphson_mean_anomaly(E, e, M):
	return E - ((M - (E - e * sin(E))) / (e * cos(E) - 1))


def newton_raphson_mean_anomaly(M, e, tol=1e-10):
	dif = float(inf)
	E, iterations = M, 1
	while dif > tol:
		E = _newton_raphson_mean_anomaly(E, e, M)
		_M = E - e * sin(E)
		dif = abs(M - _M)
		iterations += 1

	return E


def _newton_raphson_delta_eccentric_anomaly(E, r2, r2_dot, a, n, tau_i, a1):
	f = E - (1- np.linalg.norm(r2)/a) * sin(E) + a1 * (1-cos(E)) - n*tau_i
	_f = 1 - (1-np.linalg.norm(r2)/a) * cos(E) + a1 * sin(E)
	return E - (f/_f)


def newton_raphson_delta_eccentric_anomaly(tau_i, r2, r2_dot, tol=1e-14):
	dif = float(inf)
	a, n = od.OrbitalElements.calculate_semimajor_axis(r2, r2_dot), od.OrbitalElements.calculate_n(r=r2, r_dot=r2_dot)
	e = od.OrbitalElements.calculate_eccentricity(od.OrbitalElements.calculate_angular_momentum(r2, r2_dot), a)
	a1 = np.dot(r2, r2_dot) / (n * a**2)
	val = np.sign(a1 * cos(n*tau_i - a1) + (1 - np.linalg.norm(r2) / a) * sin(n * tau_i - a1))
	delta_E, iterations = n*tau_i + 0.85*e*val - a1, 1
	while dif > tol:
		_delta_E = _newton_raphson_delta_eccentric_anomaly(delta_E, r2, r2_dot, a, n, tau_i, a1)
		dif = abs(delta_E - _delta_E)
		iterations += 1
		delta_E = _delta_E
	return delta_E


if __name__ == '__main__':
	_test()
