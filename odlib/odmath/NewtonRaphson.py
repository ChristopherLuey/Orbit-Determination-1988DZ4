from math import sin, cos, inf


def test():
	M, e, tol, dif = 0.42, 0.8, 1e-10, float(inf)
	E, iterations = M, 1
	while dif > tol:
		E = _newton_raphson(E, e, M)
		_M = E - e * sin(E)
		dif = abs(M - _M)
		iterations+=1

	print("E =", E)
	print("Convergence Parameter =", tol)
	print("Iterations =", iterations)


def _newton_raphson(E, e, M):
	return E - ((M - (E - e * sin(E))) / (e * cos(E) - 1))


def newton_raphson(M, e, tol):
	dif = float(inf)
	E, iterations = M, 1
	while dif > tol:
		E = _newton_raphson(E, e, M)
		_M = E - e * sin(E)
		dif = abs(M - _M)
		iterations += 1

	return E


if __name__ == '__main__':
	test()
