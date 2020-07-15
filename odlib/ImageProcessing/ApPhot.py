import numpy as np

try:
	from astropy.io.fits import getdata
except:
	raise Exception("Cannot locate astropy")

dark = 10  # e-/pix
gain = 0.8  # e-/ADU
read = 11  # e-/pix


# ________________________
# Still a work in progress I'm having trouble with the signal being negative


def solve():
	return


def main():
	_f = "aptest.fit"
	x, y = 490.0, 293.0
	radius = 5.0
	inner_annulus, outer_annulus = 8.0, 13.0

	# _f = input("Enter a file path: ")
	# try:
	# 	x, y = list(map(float, input("X, Y: ")Newton-Raphson.split(",")))
	# 	radius = float(input("Enter a radius for the circular aperture: "))
	# 	inner_annulus, outer_annulus = list(map(float, input("Enter the inner and outer annulus radius: ").split(",")))
	# except: print("Invalid input")
	try:
		_image = getdata(filename=_f)
	except:
		print("Invalid file path")

	image = np.array(_image)
	x_centroid, y_centroid, uncert_x, uncert_y = findCentroid(image, int(x), int(y), int(radius))
	print(x_centroid, y_centroid)
	x_centroid, y_centroid = 489.99, 293.07

	# All border accepted
	arr = mask((x_centroid, y_centroid), radius, np.copy(image))
	adu_ap, n_ap = arr.sum()*read/gain, (arr > 0).sum()
	print(adu_ap, n_ap)

	arr = mask_an((x_centroid, y_centroid), inner_annulus, outer_annulus, np.copy(image))
	adu_an, n_an = arr.sum()*read/gain, (arr > 0).sum()
	print(adu_an, n_an)

	s = adu_ap - adu_an * (n_ap / n_an)
	print(s)


	print("Test File Input: ", _f)
	print("Object Near: ({}, {})".format(x, y))
	print("Aperture Radius:", radius)
	print("Annulus Inner Radius:", inner_annulus)
	print("Annulus Outer Radius:", outer_annulus)
	print("\nAll Border Pixels Accepted:")
	print("Signal = {}; SNR = {}".format(s, 1))


	return


def mask(i, radius, arr):
	a, b = i
	row, col = arr.shape
	y, x = np.ogrid[-a:row - a, -b:col - b]
	mask = x ** 2 + y ** 2 > radius ** 2
	arr[mask] = 0
	return arr


def mask_an(i, radius_small, radius_big, arr):
	a, b = i
	row, col = arr.shape
	y, x = np.ogrid[-a:row - a, -b:col - b]
	mask = x ** 2 + y **2  > radius_big ** 2
	arr[mask] = 0
	mask = x **2 + y **2 <= radius_small ** 2
	arr[mask] = 0

	return arr


def findCentroid(image, target_x, target_y, radius):
	try:
		_image = image[target_y - radius:target_y + radius + 1, target_x - radius:target_x + radius + 1]
	except IndexError:
		print("Target x,y or radius outside bounds of image")
		return 0, 0, 0, 0

	x_centroid = np.average(centroid(_image)) + target_x - 1
	y_centroid = np.average(centroid(np.rot90(_image, 1, (0, 1)))) + target_y - 1
	x_mean, y_mean = x_centroid - target_x + 1, y_centroid - target_y + 1

	x = dev(_image, x_mean)
	y = dev(np.rot90(_image, 1, (0, 1)), y_mean)

	return x_centroid, y_centroid, x, y  # <- replace with x centroid, y centroid, uncertainty in x, uncertainty in y


def dev(_image, x_mean):
	c, l = 0.0264583333, []
	for row in range(len(_image)):
		l.append([])
		for col in range(len(_image)):
			l[row].append((((col - x_mean * c) / np.sum(_image)) ** 2) * _image[row][col])
		l[row] = np.sum(l[row])
	return np.average(np.sqrt(l))


def centroid(_image):
	x, i = [], list(range(len(_image)))
	for j in _image:
		x.append(np.sum(np.multiply(_image, i)) / np.sum(_image))
	return x


if __name__ == '__main__':
	main()
