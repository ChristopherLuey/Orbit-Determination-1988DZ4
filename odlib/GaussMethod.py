import numpy as np
from odlib.odmath import ScalarEquationLagrange, NewtonRaphson
import odlib as od
from math import cos, sin, inf


def calculate_fg(tau, r2, r2_dot, series=False):
	if series:
		u = 1/(np.linalg.norm(r2) ** 2)
		z = np.dot(r2, r2_dot)/ (np.linalg.norm(r2) ** 2)
		q = (np.dot(r2_dot, r2_dot)/ (np.linalg.norm(r2) ** 2)) - u
		f = 1 - 0.5*u*tau**2 + 0.5*u*z*tau**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau**4
		g = tau - (1/6)*u*tau**3 + (1/4)*u*z*tau**4
	else:
		a, n = od.OrbitalElements.calculate_semimajor_axis(r2, r2_dot), od.OrbitalElements.calculate_n(r=r2, r_dot=r2_dot)
		delta_E = od.odmath.NewtonRaphson.newton_raphson_delta_eccentric_anomaly(tau, r2, r2_dot, tol=1.e-12)
		f, g = 1-(a/np.linalg.norm(r2))*(1-cos(delta_E)), tau + (1/n)*(sin(delta_E) - delta_E)
	return f, g


def calculate_c1c3(fg):
	#fg = [f1, g1, f3, g3]
	c1 = fg[1][1] / (fg[0][0] *fg[1][1] - fg[0][1] * fg[1][0])
	c3 = -fg[0][1] / (fg[0][0] *fg[1][1] - fg[0][1] * fg[1][0])
	return np.array([c1, c3])


def calculate_scaler_ranges(R, rhohat, c):
	D0 = np.dot(rhohat[0], rhohat[1] * rhohat[2])
	D1, D2, D3 = np.array([]), np.array([]), np.array([])
	for i in range(3):
		D1 = np.append(D1, np.dot(rhohat[2], np.cross(R[i], rhohat[1])))
		D2 = np.append(D2, np.dot(rhohat[2], np.cross(rhohat[0], R[i])))
		D3 = np.append(D3, np.dot(rhohat[0], np.cross(rhohat[1], R[i])))
	rho = np.array([])
	D = np.array([D1, D2, D3])
	for i in range(3): rho = np.append(rho, (c[0]*D[i][0] - D[i][1] + c[1]*D[i][2]) / (c[0]*D0))
	return rho


def calculate_position_vector(rho, rhohat, R):
	return rho*rhohat - R


def calculate_d1d3(fg):
	#fg = [[f1, g1], [f3, g3]]
	d1 = -fg[1][0] / (fg[0][0] * fg[1][1] - fg[0][1] * fg[1][0])
	d3 = fg[0][0] / (fg[0][0] * fg[1][1] - fg[0][1] * fg[1][0])
	return np.array([d1, d3])


def calculate_r1r3(fg, r2, r2_dot):
	#fg = [[f1, g1], [f3, g3]]
	return fg[0][0] * r2 + fg[0][1] * r2_dot, fg[1][0] * r2 + fg[1][1] * r2_dot


def calculate_velocity_vector(d1d3, r1r3):
	return d1d3[0]*r1r3[0] + d1d3[1]*r1r3[1]


def test():
	tau1 = -0.32618569435308475
	tau3 = 0.050840808143482484
	r2 = np.array([0.26640998194891174, -1.382856212643199, - 0.505199925482389])
	r2_dot = np.array([0.8439832722802604, -0.39937767878456487, 0.14200790188593015])
	rhohat = od.convert_RA_dec_rho(od.convert_HMS_degrees(*list(map(float, "17:27:15.08".split(":")))), od.convert_DMS_degrees(*list(map(float,"-24:44:53.0".split(":")))))
	R = np.array([-3.092452663004570E-02, 9.321463833793008E-01, 4.040478235741391E-01])
	fg_arr = np.array([calculate_fg(tau1, r2, r2_dot), calculate_fg(tau3, r2, r2_dot)])
	print(fg_arr)
	# c1c3_arr = calculate_c1c3(fg_arr)
	# rho = calculate_scaler_ranges(R, rhohat, c1c3_arr)
	# _r2 = calculate_position_vector(rho, rhohat, R)
	# d1d3 = calculate_d1d3(fg_arr)
	# r1r3 = calculate_r1r3(fg_arr, r2, r2_dot)
	# _r2_dot = calculate_velocity_vector(d1d3, r1r3)


def GaussMethod(R, rhohat_arr, r2, r2_dot, tau_arr, tol=1e-12):
	dif = float(inf)
	while dif > tol:
		fg_arr = np.array([calculate_fg(tau_arr[0], r2, r2_dot), calculate_fg(tau_arr[2], r2, r2_dot)])
		c1c3_arr = calculate_c1c3(fg_arr)
		rho_arr = calculate_scaler_ranges(R, rhohat_arr, c1c3_arr)
		_r2 = calculate_position_vector(rho_arr, rhohat_arr, R)
		d1d3 = calculate_d1d3(fg_arr)
		r1r3 = calculate_r1r3(fg_arr, r2, r2_dot)
		_r2_dot = calculate_velocity_vector(d1d3, r1r3)

	return


if __name__ == '__main__':
	test()
