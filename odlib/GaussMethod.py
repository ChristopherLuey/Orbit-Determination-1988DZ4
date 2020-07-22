import numpy as np
import odlib as od
from odlib.odmath import NewtonRaphson
from math import cos, sin, inf


def calculate_fg(tau, r2, r2_dot=None, series=False):
	if r2_dot is None:
		u = 1. / (np.linalg.norm(r2) ** 3)
		f = 1. - 0.5*u*tau**2
		g = tau - (1/6)*u*tau**3
	elif series is True:
		u = 1/(np.linalg.norm(r2) ** 3)
		z = np.dot(r2, r2_dot)/ (np.linalg.norm(r2) ** 2)
		q = (np.dot(r2_dot, r2_dot)/ (np.linalg.norm(r2) ** 2)) - u
		f = 1 - 0.5*u*tau**2 + 0.5*u*z*tau**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau**4
		g = tau - (1/6)*u*tau**3 + (1/4)*u*z*tau**4
	elif series is False:
		delta_E = od.odmath.NewtonRaphson.newton_raphson_delta_eccentric_anomaly(tau, r2, r2_dot, tol=1.e-14)
		f, g = 1-(od.OrbitalElements.calculate_semimajor_axis(r2, r2_dot)/np.linalg.norm(r2))*(1-cos(delta_E)), tau + (1/od.OrbitalElements.calculate_n(r=r2, r_dot=r2_dot))*(sin(delta_E) - delta_E)
	return f, g


def calculate_c1c3(fg):
	#fg = [[f1, g1], [f3, g3]]
	c1 = fg[1][1] / (fg[0][0] * fg[1][1] - fg[0][1] * fg[1][0])
	c3 = -fg[0][1] / (fg[0][0] *fg[1][1] - fg[0][1] * fg[1][0])
	return np.array([c1, c3])


def calculate_scaler_ranges(c, D0, D):
	# D = [[D11, D12, D13],
	#      [D21, D22, D23],
	#      [D31, D32, D33]]
	# c = [c1, c3]
	rho = np.empty((0,3))
	for i in range(3): rho = np.append(rho, (c[0]*D[i][0] - D[i][1] + c[1]*D[i][2]) / ([c[0], -1, c[1]][i]*D0))
	return rho


def calculate_D(rhohat, R):
	D0 = np.dot(rhohat[0], np.cross(rhohat[1], rhohat[2]))
	D1, D2, D3 = np.array([]), np.array([]), np.array([])
	for i in range(3):
		D1 = np.append(D1, np.dot(rhohat[2], np.cross(R[i], rhohat[1])))
		D2 = np.append(D2, np.dot(rhohat[2], np.cross(rhohat[0], R[i])))
		D3 = np.append(D3, np.dot(rhohat[0], np.cross(rhohat[1], R[i])))
	D = np.array([D1, D2, D3])
	return D0, D


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


def calculate_velocity_vector(d1d3, r):
	return d1d3[0]*r[0] + d1d3[1]*r[2]


def correct_light_travel(rho_arr, t_arr):
	return t_arr - (np.array(rho_arr)/od.cAU)


def GaussMethod(R_arr, rhohat_arr, r2, tau_arr, t_arr, text, D0, D, tol=1e-14):
	dif, iterations, r2_dot, _r2_dot, _rho_arr, _r2, rho_arr = float(inf), 0, None, None, None, None, None
	while dif > tol:
		fg_arr = np.array([calculate_fg(tau_arr[0], r2, r2_dot, False), calculate_fg(tau_arr[1], r2, r2_dot, False)])
		rho_arr = calculate_scaler_ranges(calculate_c1c3(fg_arr), D0, D)
		r, d1d3 = calculate_position_vector(rho_arr, rhohat_arr, R_arr), calculate_d1d3(fg_arr)
		_r2_dot = calculate_velocity_vector(d1d3, r)
		r2, r2_dot, _t_arr = r[1], _r2_dot, correct_light_travel(rho_arr, t_arr)
		tau_arr = np.array([od.convert_day_to_gaussian(_t_arr[0] - _t_arr[1]), od.convert_day_to_gaussian(_t_arr[2] - _t_arr[1]), od.convert_day_to_gaussian(_t_arr[2] - _t_arr[0])])

		if iterations is not 0:
			dif = abs(np.linalg.norm(rho_arr[1]) - np.linalg.norm(_rho_arr[1]))
			t = "\t\tIteration: " + str(iterations).rjust(2, "0") + "\n\t\t\tΔρ2 = {} AU\n\t\t\tlight-travel time = {} sec".format(dif, 2*12*3600*rho_arr[1]/od.cAU)
			text.append(t)
		_rho_arr, _r2 = rho_arr, r2
		iterations+=1
		if iterations > 999999: break

	return r2, r2_dot, rho_arr[1]
