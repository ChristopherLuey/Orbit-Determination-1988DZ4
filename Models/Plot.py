import matplotlib.pyplot as plt
import numpy as np


def histogram_plot():
        elements = np.loadtxt("io/MonteCarlo/MonteSimulation.txt")
        elements = np.reshape(elements, (-1, 10))

        fig, ax = plt.subplots(3, 2, figsize = (15,15))
        ax[0,1].set_xlabel(r'$e$')
        ax[0,1].set_title(r'Eccentricity ($e$)')

        ax[0,0].set_xlabel(r'$a$')
        ax[0,0].set_ylabel(r'Probability Density')
        ax[0,0].set_title(r'Semimajor Axis ($a$)')

        ax[1,0].set_xlabel(r'$i$')
        ax[1,0].set_ylabel(r'Probability Density')
        ax[1,0].set_title(r'Inclination ($i$)')

        ax[1,1].set_xlabel(r'$\omega$')
        ax[1,1].set_title(r'Argument of Perihelion ($\omega$)')

        ax[2,0].set_xlabel(r'$\Omega$')
        ax[2,0].set_ylabel(r'Probability Density')
        ax[2,0].set_title(r'Longitude of Ascending Node ($\Omega$)')

        ax[2,1].set_xlabel(r'$M$')
        ax[2,1].set_title(r'Mean Anomaly ($M$)')
        fig.tight_layout(pad=3.0)



        ax[0,1].hist(np.sort(elements[:,1], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)
        ax[0,0].hist(np.sort(elements[:,0], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)
        ax[1,0].hist(np.sort(elements[:,2], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)
        ax[1,1].hist(np.sort(elements[:,3], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)
        ax[2,0].hist(np.sort(elements[:,4], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)
        ax[2,1].hist(np.sort(elements[:,6], axis=None), bins=100, range=None, density=1, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=.8, log=False, color='tab:blue', label=None, stacked=False, data=None)

        fig.savefig("Distribution.jpeg", dpi=None, quality=100, optimize=True, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)

        plt.show()


def standard_deviation():
        elements = np.loadtxt("io/MonteCarlo/MonteSimulation.txt")
        elements = np.reshape(elements, (-1, 10))
        average = np.average(elements, axis=0)
        print(np.average(np.subtract(elements[:, 0], average[0]) ** 2) ** (1 / 2))
        print(np.average(np.subtract(elements[:, 1], average[1]) ** 2) ** (1 / 2))
        print(np.average(np.subtract(elements[:, 2], average[2]) ** 2) ** (1 / 2))
        print(np.average(np.subtract(elements[:, 3], average[3]) ** 2) ** (1 / 2))
        print(np.average(np.subtract(elements[:, 4], average[4]) ** 2) ** (1 / 2))
        print(np.average(np.subtract(elements[:, 6], average[6]) ** 2) ** (1 / 2))
