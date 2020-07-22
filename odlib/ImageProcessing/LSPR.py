import numpy as np
from odlib.__init__ import convert_DMS_degrees, convert_degrees_DMS, convert_HMS_degrees, convert_degrees_HMS

def LSPR(f):
    # Parse file
    file = open(f, 'r')
    star_pos = _parse(file)
    n = len(star_pos)

    # Create numpy arrays
    x, y, a, d = np.zeros((n)), np.zeros((n)), np.zeros((n)), np.zeros((n))
    for i, l in enumerate(star_pos):
        x[i], y[i], a[i], d[i] = l[0], l[1], l[2], l[3]

    # Solve for plate constants
    matrix = np.linalg.inv(np.array([[n, np.sum(x), np.sum(y)], [np.sum(x), np.sum(np.square(x)), np.sum(np.multiply(x, y))], [np.sum(y), np.sum(np.multiply(x, y)), np.sum(np.square(y))]]))
    c_a, c_d = np.array([[np.sum(a)], [np.dot(a, x)], [np.dot(a, y)]]), np.array([[np.sum(d)], [np.dot(d, x)], [np.dot(d, y)]])
    plate_a, plate_d = np.dot(matrix, c_a), np.dot(matrix, c_d)

    # Gather x,y io of object
    x_obj = float(input("x of centroid of unknown object: "))
    y_obj = float(input("y of centroid of unknown object: "))

    # Solve for a, d of object
    a_obj, d_obj = plate_a[0][0] + plate_a[1][0]*x_obj + plate_a[2][0]* y_obj, plate_d[0][0] + plate_d[1][0]*x_obj + plate_d[2][0]* y_obj

    # Solve for uncertainty
    chi_a, chi_d = 0, 0
    for i in range(1, n):
        chi_a+=((a[i] - plate_a[0][0] - plate_a[1][0]*x[i] - plate_a[2][0]*y[i]) ** 2)
        chi_d+=((d[i] - plate_d[0][0] - plate_d[1][0]*x[i] - plate_d[2][0]*y[i]) ** 2)

    uncert_a, uncert_d = (chi_a/(n-3))**(1/2), (chi_d/(n-3)) ** (1/2)
    a_obj, d_obj = convert_degrees_HMS(a_obj), convert_degrees_DMS(d_obj)

    # Print outputs
    # test(f, x_obj, y_obj, uncert_a, uncert_d, a_obj, d_obj, plate_a, plate_d)

    # [b1, b2, a11, a12, a21, a22], [RA of obj, dec of obj], [RA uncertainty, dec uncertainty]
    return [plate_a[0][0], plate_d[0][0], plate_a[1][0], plate_a[2][0], plate_d[1][0], plate_d[2][0]], [a_obj, d_obj], [uncert_a, uncert_d]


def _parse(file):
    star_pos = []
    for i, line in enumerate(file, 0):
        star_pos.append(line.strip("\n").split(" "))
        for j in range(len(star_pos[i])):
            if j <= 1:
                star_pos[i][j] = float(star_pos[i][j])
            elif j==3:
                a = star_pos[i][j].split(":")
                star_pos[i][j] = convert_DMS_degrees(float(a[0]), float(a[1]), float(a[2]))
            elif j==2:
                a = star_pos[i][j].split(":")
                star_pos[i][j] = convert_HMS_degrees(float(a[0]), float(a[1]), float(a[2]))
    return star_pos


def test(f, x_obj, y_obj, uncert_a, uncert_d, a_obj, d_obj, plate_a, plate_d):
    # Print outputs
    print("\n\n\n\ntest io file:", f)
    print("test position (x,y):", x_obj, y_obj)
    print("\n\n***************\nPlate Constants\n***************")
    txt = "b1: {b1} deg\nb2: {b2} deg\na11: {a11} deg/pix\na12: {a12} deg/pix\na21: {a21} deg/pix\na22: {a22} deg/pix".format(b1=plate_a[0][0], b2=plate_d[0][0], a11=plate_a[1][0], a12=plate_a[2][0], a21=plate_d[1][0], a22=plate_d[2][0])
    print(txt)
    print("\n\n***********\nUncertainty\n***********")
    print("RA: {} arcseconds\nDec: {} arcseconds".format(uncert_a*3600, uncert_d*3600))
    print("\n\n*********************************\nAstrometry for\n(x,y): ({},{})".format(x_obj, y_obj))
    print("*********************************")
    if d_obj[0] < 0:print("RA: {}:{}:{}\nDec: {}{}:{}:{}".format(a_obj[0], a_obj[1], a_obj[2], "-", abs(d_obj[0]),d_obj[1], d_obj[2]))
    else:print("RA: {}:{}:{}\nDec: {}{}:{}:{}".format(a_obj[0], a_obj[1], a_obj[2], "+", abs(d_obj[0]),d_obj[1], d_obj[2]))


if __name__ == '__main__':
    test()