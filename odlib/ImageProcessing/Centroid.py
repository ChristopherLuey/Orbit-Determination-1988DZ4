# Programming Assignment #3 / SSP
# Christopher Luey
# 6/30

import numpy as np
from astropy.io.fits import getdata


def findCentroid(fits_file, target_x, target_y, radius):
    try:
        image = getdata(filename=fits_file)
        _image = image[target_y-radius:target_y+radius+1, target_x-radius:target_x+radius+1]
    except FileNotFoundError:
        print("Could not find file:", fits_file)
        return 0, 0, 0, 0
    except IndexError:
        print("Target x,y or radius outside bounds of image")
        return 0, 0, 0, 0

    x_centroid = np.average(centroid(_image)) + target_x - 1
    y_centroid = np.average(centroid(np.rot90(_image, 1, (0,1)))) + target_y - 1
    x_mean, y_mean = x_centroid - target_x + 1, y_centroid - target_y + 1

    x = dev(_image, x_mean)
    y = dev(np.rot90(_image, 1, (0,1)), y_mean)

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


# centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("../sampleimage.fits", 351, 154, 1)
# if abs(centroid_x - 350.9958) < 1e-3 and abs(centroid_y - 153.9955) < 1e-3:
#     print("centroid calculation CORRECT")
# else:
#     print(f"centroid calculation INCORRECT, expected (350.9958, 153.9955), got ({centroid_x}, {centroid_y})")
# if abs(uncert_x - 0.005254018) < 1e-6 and abs(uncert_y - 0.005249733) < 1e-6:
#     print("uncertainty calculation CORRECT")
# else:
#     print(f"uncertainty calculation INCORRECT, expected (0.005254018, 0.005249733), got ({uncert_x}, {uncert_y})")
